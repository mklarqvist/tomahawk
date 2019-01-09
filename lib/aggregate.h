/*
Copyright (C) 2016-present Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/
#include <getopt.h>

#include "utility.h"
#include "two_reader.h"
#include "sort_progress.h"

namespace tomahawk {

struct offset_tuple {
    offset_tuple() : range(0), min(std::numeric_limits<uint32_t>::max()), max(0){}
    uint64_t range;
    uint32_t min, max;
};

/**<
 *
 */
// Todo:
struct twk1_aggregate {
public:
	twk1_aggregate() : n(0), x(0), y(0), bpx(0), bpy(0), n_original(0), range(0), data(nullptr){}
	twk1_aggregate(const uint32_t x, const uint32_t y) : n(x*y), x(x), y(y), bpx(0), bpy(0), n_original(0), range(0), data(new double[n]){}
	~twk1_aggregate(){ delete[] data; }

	friend std::ostream& operator<<(std::ostream& stream, const twk1_aggregate& agg){
		stream.write(TOMAHAWK_AGGREGATE_MAGIC_HEADER.data(), TOMAHAWK_AGGREGATE_MAGIC_HEADER_LENGTH);

		SerializePrimitive(agg.n, stream);
		SerializePrimitive(agg.x, stream);
		SerializePrimitive(agg.y, stream);
		SerializePrimitive(agg.bpx, stream);
		SerializePrimitive(agg.bpy, stream);
		SerializePrimitive(agg.n_original, stream);
		SerializePrimitive(agg.range, stream);
		SerializeString(agg.filename, stream);

		// Write rid tuples.
		uint32_t n_rid = agg.rid_offsets.size();
		SerializePrimitive(n_rid, stream);
		for(int i = 0; i < agg.rid_offsets.size(); ++i){
			uint32_t min = agg.rid_offsets[i].min == std::numeric_limits<uint32_t>::max() ? 0 : agg.rid_offsets[i].min;
			uint32_t max = agg.rid_offsets[i].max < min ? 0 : agg.rid_offsets[i].min;
			SerializePrimitive(min, stream);
			SerializePrimitive(max, stream);
			SerializePrimitive(agg.rid_offsets[i].range, stream);
		}

		ZSTDCodec zcodec;
		twk_buffer_t ibuf, obuf;
		for(int i = 0; i < agg.n; ++i) ibuf += agg.data[i];
		zcodec.Compress(ibuf, obuf, 6);

		// Write data.
		uint32_t obuf_size = obuf.size();
		SerializePrimitive(obuf_size, stream);
		stream.write(obuf.data(), obuf.size());
		stream.write(TOMAHAWK_TWOAGG_EOF.data(), TOMAHAWK_TWOAGG_EOF_LENGTH);

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, twk1_aggregate& agg){
		char magic[TOMAHAWK_AGGREGATE_MAGIC_HEADER_LENGTH];
		stream.read(&magic[0], TOMAHAWK_AGGREGATE_MAGIC_HEADER_LENGTH);

		DeserializePrimitive(agg.n, stream);
		DeserializePrimitive(agg.x, stream);
		DeserializePrimitive(agg.y, stream);
		DeserializePrimitive(agg.bpx, stream);
		DeserializePrimitive(agg.bpy, stream);
		DeserializePrimitive(agg.n_original, stream);
		DeserializePrimitive(agg.range, stream);
		DeserializeString(agg.filename, stream);

		// Write rid tuples.
		uint32_t n_rid = 0;
		DeserializePrimitive(n_rid, stream);
		agg.rid_offsets.clear();
		agg.rid_offsets.resize(n_rid);
		for(int i = 0; i < agg.rid_offsets.size(); ++i){
			DeserializePrimitive(agg.rid_offsets[i].min, stream);
			DeserializePrimitive(agg.rid_offsets[i].max, stream);
			DeserializePrimitive(agg.rid_offsets[i].range, stream);
		}

		ZSTDCodec zcodec;
		twk_buffer_t ibuf, obuf;
		for(int i = 0; i < agg.n; ++i) ibuf += agg.data[i];
		zcodec.Decompress(obuf, ibuf);

		// Write data.
		delete[] agg.data; agg.data = nullptr;
		agg.data = new double[agg.n];
		uint32_t obuf_size = 0;
		DeserializePrimitive(obuf_size, stream);
		stream.read(obuf.data(), obuf_size);

		char eof[TOMAHAWK_TWOAGG_EOF_LENGTH];
		stream.read(&eof[0], TOMAHAWK_TWOAGG_EOF_LENGTH);

		return(stream);
	}

	bool Open(std::string input){
		if(input.size() == 0) return false;

		std::ifstream in(input, std::ios::in | std::ios::binary | std::ios::ate);
		if(in.good() == false) return false;
		uint64_t fsize = in.tellg();
		in.seekg(0);

		char magic[TOMAHAWK_AGGREGATE_MAGIC_HEADER_LENGTH];
		in.read(&magic[0], TOMAHAWK_AGGREGATE_MAGIC_HEADER_LENGTH);
		if(in.good() == false) return false;
		if(strncmp(magic, TOMAHAWK_AGGREGATE_MAGIC_HEADER.data(), TOMAHAWK_AGGREGATE_MAGIC_HEADER_LENGTH) != 0) return false;

		DeserializePrimitive(n, in);
		DeserializePrimitive(x, in);
		DeserializePrimitive(y, in);
		DeserializePrimitive(bpx, in);
		DeserializePrimitive(bpy, in);
		DeserializePrimitive(n_original, in);
		DeserializePrimitive(range, in);
		DeserializeString(filename, in);
		if(in.good() == false) return false;

		// Write rid tuples.
		uint32_t n_rid = 0;
		DeserializePrimitive(n_rid, in);
		rid_offsets.clear();
		rid_offsets.resize(n_rid);
		for(int i = 0; i < rid_offsets.size(); ++i){
			DeserializePrimitive(rid_offsets[i].min, in);
			DeserializePrimitive(rid_offsets[i].max, in);
			DeserializePrimitive(rid_offsets[i].range, in);
		}
		if(in.good() == false) return false;

		ZSTDCodec zcodec;
		twk_buffer_t ibuf, obuf;
		for(int i = 0; i < n; ++i) ibuf += data[i];
		zcodec.Decompress(obuf, ibuf);

		// Write data.
		delete[] data; data = nullptr;
		data = new double[n];
		uint32_t obuf_size = 0;
		DeserializePrimitive(obuf_size, in);
		in.read(obuf.data(), obuf_size);
		assert(in.tellg() == fsize - TOMAHAWK_TWOAGG_EOF_LENGTH);

		//in.seekg(fsize - TOMAHAWK_TWOAGG_EOF_LENGTH);
		if(in.good() == false) return false;
		char eof[TOMAHAWK_TWOAGG_EOF_LENGTH];
		in.read(&eof[0], TOMAHAWK_TWOAGG_EOF_LENGTH);
		if(strncmp(eof, TOMAHAWK_TWOAGG_EOF.data(), TOMAHAWK_TWOAGG_EOF_LENGTH) != 0) return false;


		return false;
	}

public:
	// magic header
	uint32_t n, x, y, bpx, bpy, n_original;
	uint64_t range;
	std::string filename; // input filename
	std::vector<offset_tuple> rid_offsets; // mat offsets
	double* data;
	// EOF
};

struct twk_agg_slave {
public:
	struct range_helper {
		range_helper() : set(false), min(std::numeric_limits<uint32_t>::max()), max(0){}

		bool set;
		uint32_t min, max;
	};

public:
	twk_agg_slave() : f(0), t(0), xrange(0), yrange(0), it(nullptr),
	    thread(nullptr), progress(nullptr),
	    aggregator(&twk_sstats::AddR2), reductor(&twk_sstats::GetMean)
	{}
	~twk_agg_slave(){ delete it; delete thread; }

	/**<
	 *
	 * @param rdr Reference instance of a two reader.
	 * @return    Returns a pointer to the spawned thread instance if successful or a nullptr otherwise.
	 */
	std::thread* StartFindRanges(two_reader& rdr){
		if(f > t) return nullptr;
		if(t - f == 0) return nullptr;

		if(stream.good()) stream.close();

		stream = std::ifstream(filename,std::ios::binary | std::ios::in);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to open \"" << filename << "\"..." << std::endl;
			return nullptr;
		}

		stream.seekg(rdr.index.ent[f].foff);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to seek to position " << rdr.index.ent[f].foff << " in \""  << filename << "\"..." << std::endl;
			return nullptr;
		}

		delete it; it = nullptr;
		it  = new twk1_two_iterator;
		it->stream = &stream;

		// Vector of contigs.
		contig_avail.resize(rdr.hdr.GetNumberContigs());

		// New thread.
		delete thread; thread = nullptr;
		thread = new std::thread(&twk_agg_slave::FindRangesUnsorted, this);

		return(thread);
	}

	/**<
	 *
	 * @param rdr Reference instance of a two reader.
	 * @param agg Aggregation function pointer.
	 * @param red Reduction function pointer.
	 * @param x   Number of bins in X-dimension.
	 * @param y   Number of bins in Y-dimension.
	 * @param xr  Range of X.
	 * @param yr  Range of Y.
	 * @return    Returns a pointer to the spawned thread instance if successful or a nullptr otherwise.
	 */
	std::thread* StartBuildMatrix(two_reader& rdr,
			twk_sstats::aggfunc agg,
			twk_sstats::redfunc red,
			uint32_t x, uint32_t y,
			uint32_t xr, uint32_t yr)
	{
		if(f > t) return nullptr;
		if(t - f == 0) return nullptr;

		assert(agg != nullptr);
		assert(red != nullptr);
		aggregator = agg;
		reductor = red;
		xrange = xr;
		yrange = yr;
		mat = std::vector< std::vector<twk_sstats> >(x, std::vector<twk_sstats>(y));

		if(stream.good()) stream.close();

		stream = std::ifstream(filename,std::ios::binary | std::ios::in);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to open \"" << filename << "\"..." << std::endl;
			return nullptr;
		}

		stream.seekg(rdr.index.ent[f].foff);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to seek to position " << rdr.index.ent[f].foff << " in \""  << filename << "\"..." << std::endl;
			return nullptr;
		}

		delete it; it = nullptr;
		it  = new twk1_two_iterator;
		it->stream = &stream;

		// New thread.
		delete thread; thread = nullptr;
		thread = new std::thread(&twk_agg_slave::BuildMatrix, this);

		return(thread);
	}

	/**<
	 * Investigate the presence of blocks mapping to unique chromosomes. This
	 * information is needed to set the (X,Y)-coordinate system for the aggregation
	 * space.
	 * @param rdr Input reference reader.
	 * @return    Returns TRUE upon success or FALSE otherwise.
	 */
	bool FindRangesUnsorted(){
		uint32_t tot = 0;
		for(int i = f; i < t; ++i){
			assert(it->NextBlock());
			// Todo: bugfix this case
			if(it->blk.n < 5) continue;

			tot += it->GetBlock().n;
			//std::cerr << "here" << std::endl;

			for(int j = 0; j < it->blk.n; ++j){
				//std::cerr << j << "/" << it->blk.n << "\t" << it->blk[j].ridA << "/" << contig_avail.size() << std::endl;
				assert(it->blk[j].ridA < contig_avail.size());
				assert(it->blk[j].ridB < contig_avail.size());

				contig_avail[it->blk[j].ridA].set = true;
				contig_avail[it->blk[j].ridB].set = true;
				contig_avail[it->blk[j].ridA].min = std::min(it->blk[j].Apos, contig_avail[it->blk[j].ridA].min);
				contig_avail[it->blk[j].ridA].max = std::max(it->blk[j].Apos, contig_avail[it->blk[j].ridA].max);
				contig_avail[it->blk[j].ridB].min = std::min(it->blk[j].Bpos, contig_avail[it->blk[j].ridB].min);
				contig_avail[it->blk[j].ridB].max = std::max(it->blk[j].Bpos, contig_avail[it->blk[j].ridB].max);
			}
			progress->cmps += it->GetBlock().n;
		}

		delete it; it = nullptr;
		return true;
	}

	bool BuildMatrix(){
		uint32_t tot = 0;
		for(int i = f; i < t; ++i){
			assert(it->NextBlock());
			// Todo: bugfix this case
			if(it->blk.n < 5) continue;

			tot += it->GetBlock().n;
			for(int j = 0; j < it->blk.n; ++j){
				// Invoke aggregator function.
				// Position: cumulative offset up to chromosome + left-adjusted position
				// Position: (chromosome_offset.range - chromosome_offset.max - chromoosme_offset.min) + (Apos - smallest_in_chr)
				(mat[((rid_offsets[it->blk[j].ridA].range - (rid_offsets[it->blk[j].ridA].max - rid_offsets[it->blk[j].ridA].min)) + (it->blk[j].Apos - rid_offsets[it->blk[j].ridA].min))/xrange][((rid_offsets[it->blk[j].ridB].range - (rid_offsets[it->blk[j].ridB].max - rid_offsets[it->blk[j].ridB].min)) + (it->blk[j].Bpos - rid_offsets[it->blk[j].ridB].min))/yrange].*aggregator)(&it->blk[j]);
			}
			progress->cmps += it->GetBlock().n;
		}

		delete it; it = nullptr;
		return true;
	}

	/**<
	 * Reduction helper for matrices. Adds the data from the matrix in the passed
	 * slave instance to the current instance.
	 * @param other Reference to other slave instance.
	 */
	void AddMatrix(const twk_agg_slave& other){
		for(int i = 0; i < mat.size(); ++i){
			for(int j = 0; j < mat[i].size(); ++j)
				mat[i][j] += other.mat[i][j];
		}
	}

	/**<
	 * Simple print subroutine for the internal matrix.
	 * @param stream     Target output stream to print to.
	 * @param min_cutoff Minimum number of elements require to output a value. If the observed frequency is smaller than this value then we emit 0.
	 */
	void PrintMatrix(std::ostream& stream, uint32_t min_cutoff = 5) const{
		for(int i = 0; i < mat.size(); ++i){
			stream << (mat[i][0].*reductor)(min_cutoff);
			for(int j = 1; j < mat[i].size(); ++j)
				stream << '\t' << (mat[i][j].*reductor)(min_cutoff);
			stream << '\n';
		}
		stream.flush();
	}

	twk1_aggregate& Overload(twk1_aggregate& agg, uint32_t min_cutoff = 5) const{
		for(int i = 0; i < mat.size(); ++i){
			agg.data[i*mat.size() + 0] = (mat[i][0].*reductor)(min_cutoff);
			for(int j = 1; j < mat[i].size(); ++j){
				agg.data[i*mat.size() + j] = (mat[i][j].*reductor)(min_cutoff);
				//stream << '\t' << (mat[i][j].*reductor)(min_cutoff);
			}
			//stream << '\n';
		}
		//stream.flush();
		return(agg);
	}

public:
	uint32_t f, t; // (from,to)-tuple
	uint32_t xrange, yrange; // (x,y)-tuple
	std::ifstream stream;
	twk1_two_iterator* it;
	std::thread* thread;
	twk_sort_progress* progress;
	twk_sstats::aggfunc aggregator; // aggregator function
	twk_sstats::redfunc reductor; // reductor function
	std::string filename; // input filename
	std::vector<range_helper> contig_avail;
	std::vector<offset_tuple> rid_offsets; // mat offsets
	std::vector< std::vector<twk_sstats> > mat; // Output matrix
};

}

void aggregate_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Aggregate TWO data into a rasterized matrix of size [x,y] for\n"
    "        plotting.\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " aggregate [options] <in.two>\n\n"

	"Options:\n"
	"  -i   FILE   input TWO file (required)\n"
	"  -x,y INT    number of X/Y-axis bins (default: 1000)\n"
	"  -f   STRING aggregation function: can be one of (r2,r,d,dprime,dp,p,hets,alts,het,alt)(required)\n"
	"  -r   STRING reduction function: can be one of (mean,count,n,min,max,sd)(required)\n"
	"  -I   STRING filter interval <contig>:pos-pos (TWK/TWO) or linked interval <contig>:pos-pos,<contig>:pos-pos\n"
	"  -c   INT    min cut-off value used in reduction function: value < c will be set to 0 (default: 5)\n"
	"  -t   INT    number of parallel threads: each thread will use " << sizeof(tomahawk::twk_sstats) << "(x*y) bytes\n" << std::endl;
}

int aggregate(int argc, char** argv){
	if(argc < 3){
		aggregate_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{"interval",    optional_argument, 0, 'I' },
		{"xbins",    optional_argument, 0, 'x' },
		{"ynins",    optional_argument, 0, 'y' },
		{"min-cutoff",    optional_argument, 0, 'c' },
		{"aggregate-function",    required_argument, 0, 'f' },
		{"reduce-function",    required_argument, 0, 'r' },
		{"threads",    optional_argument, 0, 't' },
		{0,0,0,0}
	};

	tomahawk::twk_two_settings settings;
	int32_t x_bins = 1000, y_bins = 1000;
	std::string aggregate_func_name, reduce_func_name;
	int32_t min_cutoff = 5;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:I:x:y:f:r:c:t:", long_options, &long_index)) != -1){
		hits += 2;
		switch (c){
		case ':':   /* missing option argument */
			fprintf(stderr, "%s: option `-%c' requires an argument\n",
					argv[0], optopt);
			break;

		case '?':
		default:
			fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
					argv[0], optopt);
			break;

		case 'i':
			settings.in = std::string(optarg);
			break;
		case 'I': settings.ivals.push_back(std::string(optarg)); break;
		case 'x': x_bins = std::atoi(optarg); break;
		case 'y': y_bins = std::atoi(optarg); break;
		case 'f': aggregate_func_name = std::string(optarg); break;
		case 'r': reduce_func_name = std::string(optarg); break;
		case 'c': min_cutoff = std::atoi(optarg); break;
		case 't': settings.n_threads = std::atoi(optarg); break;
		}
	}

	tomahawk::twk_sstats::aggfunc f = &tomahawk::twk_sstats::AddR2;
	tomahawk::twk_sstats::redfunc r = &tomahawk::twk_sstats::GetMean;

	if(aggregate_func_name.size() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No aggregation function (-f) provided..." << std::endl;
		return(1);
	}

	if(reduce_func_name.size() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No reduce function (-r) provided..." << std::endl;
		return(1);
	}

	if(min_cutoff < 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a min-cutoff (-c) < 0..." << std::endl;
		return(1);
	}

	if(settings.n_threads <= 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have <= 0 threads (-t)..." << std::endl;
		return(1);
	}

	// Transform string of aggregation function name into lower then try to map
	// name to existing function names.
	std::transform(aggregate_func_name.begin(), aggregate_func_name.end(), aggregate_func_name.begin(), ::tolower);
	if(aggregate_func_name == "r2")         { f = &tomahawk::twk_sstats::AddR2;   }
	else if(aggregate_func_name == "r")     { f = &tomahawk::twk_sstats::AddR;    }
	else if(aggregate_func_name == "d")     { f = &tomahawk::twk_sstats::AddD;    }
	else if(aggregate_func_name == "dprime"){ f = &tomahawk::twk_sstats::AddDprime; }
	else if(aggregate_func_name == "dp")    { f = &tomahawk::twk_sstats::AddDprime; }
	else if(aggregate_func_name == "p")     { f = &tomahawk::twk_sstats::AddP;    }
	else if(aggregate_func_name == "hets")  { f = &tomahawk::twk_sstats::AddHets; }
	else if(aggregate_func_name == "alts")  { f = &tomahawk::twk_sstats::AddAlts; }
	else if(aggregate_func_name == "het")   { f = &tomahawk::twk_sstats::AddHets; }
	else if(aggregate_func_name == "alt")   { f = &tomahawk::twk_sstats::AddAlts; }
	else {
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Unknown aggregation function \"" << aggregate_func_name << "\"..." << std::endl;
		return(1);
	}

	// Transform string of reduction function name into lower then try to map
    // name to existing function names.
	std::transform(reduce_func_name.begin(), reduce_func_name.end(), reduce_func_name.begin(), ::tolower);
	if(reduce_func_name == "mean")          { r = &tomahawk::twk_sstats::GetMean;  }
		else if(reduce_func_name == "max")  { r = &tomahawk::twk_sstats::GetMax;   }
		else if(reduce_func_name == "min")  { r = &tomahawk::twk_sstats::GetMin;   }
		else if(reduce_func_name == "count"){ r = &tomahawk::twk_sstats::GetTotal; }
		else if(reduce_func_name == "n")    { r = &tomahawk::twk_sstats::GetTotal; }
		else if(reduce_func_name == "sd")   { r = &tomahawk::twk_sstats::GetStandardDeviation; }
		else {
			std::cerr << tomahawk::utility::timestamp("ERROR") << "Unknown reduce function \"" << reduce_func_name << "\"..." << std::endl;
			return(1);
		}

	if(x_bins < 5){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Number of x-bins cannot be < 5!" << std::endl;
		return(1);
	}

	if(y_bins < 5){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Number of y-bins cannot be < 5!" << std::endl;
		return(1);
	}

	if(settings.in.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	// New instance of reader.
	tomahawk::two_reader oreader;

	// Open file handle.
	if(oreader.Open(settings.in) == false) return 1;


	// Build intervals data structures if any are available.
	/*if(settings.intervals.Build(settings.ivals,
	                            oreader.hdr.GetNumberContigs(),
	                            oreader.index,
	                            oreader.hdr) == false)
	{
		return 1;
	}*/

	// Construct filters.
	//settings.filter.Build();

	// Print messages
	tomahawk::ProgramMessage();
	std::cerr << tomahawk::utility::timestamp("LOG") << "Calling aggregate..." << std::endl;


	// Algorithmic overview.
	// Step 1: Find maximum and minimum X and Y values.
	// Step 2: Partition into (maxY-minY)/#bins and (maxX-minX)/#bins buckets.
	//      a: Range can be either dynamic min and max given the data, or;
	//      b: Interval (from,to)-tuple.
	// Step 3: Iterate over data and update summary statistics in buckets.
	// Step 4: Output data.
	std::cerr << tomahawk::utility::timestamp("LOG") << "Performing 2-pass over data..." << std::endl;

	// Step 1: iterate over index entries and find what contigs are used
	std::cerr << tomahawk::utility::timestamp("LOG") << "===== First pass (peeking at landscape) =====" << std::endl;
	// Distrubution.
	uint64_t b_unc = 0, n_recs = 0;
	std::cerr << tomahawk::utility::timestamp("LOG") << "Blocks: " << tomahawk::utility::ToPrettyString(oreader.index.n) << std::endl;
	for(int i = 0; i < oreader.index.n; ++i){
		b_unc  += oreader.index.ent[i].b_unc;
		n_recs += oreader.index.ent[i].n;
	}
	std::cerr << tomahawk::utility::timestamp("LOG") << "Uncompressed size: " << tomahawk::utility::ToPrettyDiskString(b_unc) << std::endl;

	if(b_unc == 0){
		std::cerr << tomahawk::utility::timestamp("LOG") << "Cannot aggregate empty file..." << std::endl;
		return false;
	}

	if(oreader.index.n < settings.n_threads) settings.n_threads = oreader.index.n;
	uint64_t b_unc_thread = b_unc / settings.n_threads;
	std::cerr << tomahawk::utility::timestamp("LOG","THREAD") << "Data/thread: " << tomahawk::utility::ToPrettyDiskString(b_unc_thread) << std::endl;

	std::vector< std::pair<uint32_t,uint32_t> > ranges;
	uint64_t fR = 0, tR = 0, b_unc_tot = 0;
	for(int i = 0; i < oreader.index.n; ++i){
		if(b_unc_tot >= b_unc_thread){
			ranges.push_back(std::pair<uint32_t,uint32_t>(fR, tR));
			b_unc_tot = 0;
			fR = tR;
		}
		b_unc_tot += oreader.index.ent[i].b_unc;
		++tR;
	}
	if(fR != tR){
		ranges.push_back(std::pair<uint32_t,uint32_t>(fR, tR));
		b_unc_tot = 0;
		fR = tR;
	}
	assert(ranges.back().second == oreader.index.n);
	assert(ranges.size() <= settings.n_threads);

	tomahawk::twk_sort_progress progress_sort;
	progress_sort.n_cmps = n_recs;
	std::thread* psthread = progress_sort.Start();

	tomahawk::twk_agg_slave* slaves = new tomahawk::twk_agg_slave[settings.n_threads];
	uint32_t range_thread = oreader.index.n / settings.n_threads;
	for(int i = 0; i < settings.n_threads; ++i){
		slaves[i].f = ranges[i].first;
		slaves[i].t = ranges[i].second;
		slaves[i].filename = settings.in;
		slaves[i].progress = &progress_sort;
		//std::cerr << "thread-" << i << " " << slaves[i].f << "-" << slaves[i].t << std::endl;
	}

	for(int i = 0; i < settings.n_threads; ++i){
		if(slaves[i].StartFindRanges(oreader) == nullptr){
			std::cerr << tomahawk::utility::timestamp("ERROR","THREAD") << "Failed to spawn slave" << std::endl;
			return false;
		}
	}
	for(int i = 0; i < settings.n_threads; ++i) slaves[i].thread->join();
	progress_sort.is_ticking = false;
	progress_sort.PrintFinal();

	// Reduce.
	for(int i = 1; i < settings.n_threads; ++i){
		for(int j = 0; j < slaves[0].contig_avail.size(); ++j){
			slaves[0].contig_avail[j].set = std::max(slaves[0].contig_avail[j].set, slaves[i].contig_avail[j].set);
			slaves[0].contig_avail[j].min = std::min(slaves[0].contig_avail[j].min, slaves[i].contig_avail[j].min);
			slaves[0].contig_avail[j].max = std::max(slaves[0].contig_avail[j].max, slaves[i].contig_avail[j].max);
		}
	}

	// Reduce.
	uint32_t n_chrom_set = slaves[0].contig_avail[0].set;
	for(int i = 0; i < slaves[0].contig_avail.size(); ++i){
		n_chrom_set += slaves[0].contig_avail[i].set;
	}
	//std::cerr << "chrom set=" << n_chrom_set << std::endl;

	// Step 2: Determine boundaries given the contigs that were set.
	//         Calculate the landscape ranges (X and Y dimensions).
	// Approach 2: Dropping regions with no data.
	uint64_t range = 0;
	std::vector<tomahawk::offset_tuple> rid_offsets(slaves[0].contig_avail.size());

	/**<
	 * If there is only chromosome set then restrict the (X,Y) landscape to the
	 * available data range. This is in contrast to cases where N > 1, where we
	 * consider the entire genomic range of the affected chromosomes irrespective
	 * of how much range is actually used.
	 */
	if(n_chrom_set == 1){
		if(slaves[0].contig_avail[0].set){
			rid_offsets[0].range = slaves[0].contig_avail[0].max - slaves[0].contig_avail[0].min + 1;
			range += slaves[0].contig_avail[0].max - slaves[0].contig_avail[0].min + 1;
		} else
			rid_offsets[0].range = 0;

		rid_offsets[0].min = slaves[0].contig_avail[0].min;
		rid_offsets[0].max = slaves[0].contig_avail[0].max;

		for(int i = 1; i < slaves[0].contig_avail.size(); ++i){
			if(slaves[0].contig_avail[i].set){
				// Cumulative offset for the current rid equals the previous rid
				rid_offsets[i].range = rid_offsets[i - 1].range + (slaves[0].contig_avail[i].max - slaves[0].contig_avail[i].min + 1);
				range += (slaves[0].contig_avail[i].max - slaves[0].contig_avail[i].min + 1);
			} else {
				rid_offsets[i].range = rid_offsets[i - 1].range;
			}
			rid_offsets[i].min = slaves[0].contig_avail[i].min;
			rid_offsets[i].max = slaves[0].contig_avail[i].max;
		}

	}
	// If there is data from n>1 chromosomes.
	else {
		if(slaves[0].contig_avail[0].set){
			rid_offsets[0].range = oreader.hdr.contigs_[0].n_bases;
			range += oreader.hdr.contigs_[0].n_bases;
		} else
			rid_offsets[0].range = 0;

		rid_offsets[0].min = 0;
		rid_offsets[0].max = oreader.hdr.contigs_[0].n_bases;

		for(int i = 1; i < slaves[0].contig_avail.size(); ++i){
			if(slaves[0].contig_avail[i].set){
				// Cumulative offset for the current rid equals the previous rid
				rid_offsets[i].range = rid_offsets[i - 1].range + oreader.hdr.contigs_[i].n_bases;
				range += oreader.hdr.contigs_[i].n_bases;

			} else {
				rid_offsets[i].range = rid_offsets[i - 1].range;
			}
			rid_offsets[i].min = 0;
			rid_offsets[i].max = oreader.hdr.contigs_[i].n_bases;
		}
	}

	//std::cerr << "range=" << range << std::endl;
	//for(int i = 0; i < rid_offsets.size(); ++i){
	//	std::cerr << "rid=" << i << "=" << rid_offsets[i].range << " -> " << rid_offsets[i].min << "-" << rid_offsets[i].max << std::endl;
	//}

	// Step 3: Second pass over data.
	//         Prepare n-tensor for storing output data.
	//         Matrix dimensions (1,2) correspond to pixel equivalents.
	//         Tensor dimensions (3,..) correspond to summary statistics for
	//         each bin (pixel).
	uint32_t xrange = std::ceil((float)range / x_bins);
	uint32_t yrange = std::ceil((float)range / y_bins);

	tomahawk::twk1_aggregate agg(x_bins, y_bins);
	agg.bpx = xrange;
	agg.bpy = yrange;
	agg.range = range;
	agg.rid_offsets = rid_offsets;
	agg.n_original = n_recs;

	std::cerr << tomahawk::utility::timestamp("LOG") << "===== Second pass (building matrix) =====" << std::endl;
	std::cerr << tomahawk::utility::timestamp("LOG") << "Aggregating " << tomahawk::utility::ToPrettyString(n_recs) << " records..." << std::endl;
	std::cerr << tomahawk::utility::timestamp("LOG","THREAD") << "Allocating: " << tomahawk::utility::ToPrettyDiskString(sizeof(tomahawk::twk_sstats)*x_bins*y_bins*settings.n_threads) << " for matrices..." << std::endl;

	tomahawk::twk_sort_progress progress_sort_step2;
	progress_sort_step2.n_cmps = n_recs;
	psthread = progress_sort_step2.Start();
	//std::cerr << "range=" << range << " x,y = " << xrange << " bp/pixel " << " and " << yrange << " bp/pixel" << std::endl;
	for(int i = 0; i < settings.n_threads; ++i){
		slaves[i].rid_offsets = rid_offsets;
		slaves[i].progress = &progress_sort_step2;
		if(slaves[i].StartBuildMatrix(oreader, f, r, x_bins, y_bins, xrange, yrange) == nullptr){
			std::cerr << tomahawk::utility::timestamp("ERROR","THREAD") << "Failed to spawn slave" << std::endl;
			return false;
		}
	}
	for(int i = 0; i < settings.n_threads; ++i) slaves[i].thread->join();
	progress_sort_step2.is_ticking = false;
	progress_sort_step2.PrintFinal();
	for(int i = 1; i < settings.n_threads; ++i) slaves[0].AddMatrix(slaves[i]);

	// Print matrix
	//slaves[0].PrintMatrix(std::cout, min_cutoff);
	slaves[0].Overload(agg);
	std::cout << agg;

	std::cerr << tomahawk::utility::timestamp("LOG") << "Aggregated " << tomahawk::utility::ToPrettyString(n_recs) << " records in " << tomahawk::utility::ToPrettyString(x_bins*y_bins) << " bins." << std::endl;
	std::cerr << tomahawk::utility::timestamp("LOG") << "Finished." << std::endl;

	delete[] slaves;

	return(0);
}
