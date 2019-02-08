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
#include "aggregation.h"

void aggregate_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Aggregate TWO data into a rasterized matrix of size [x,y] for\n"
    "        plotting.\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " aggregate [options] <in.two>\n\n"

	"Options:\n"
	"  -i   FILE   input TWO file (required)\n"
	"  -o   FILE   output file path (default: -)\n"
	"  -O   <b,u>  b: compressed binary representation, u: uncompressed matrix (default: b)\n"
	"  -x,y INT    number of X/Y-axis bins (default: 1000)\n"
	"  -f   STRING aggregation function: can be one of (r2,r,d,dprime,dp,p,hets,alts,het,alt)(required)\n"
	"  -r   STRING reduction function: can be one of (mean,count,n,min,max,sd,total)(required)\n"
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
		{"output",       optional_argument, 0, 'o' },
		{"output-type",       optional_argument, 0, 'O' },
		{"interval",    optional_argument, 0, 'I' },
		{"xbins",       optional_argument, 0, 'x' },
		{"ynins",       optional_argument, 0, 'y' },
		{"min-cutoff",  optional_argument, 0, 'c' },
		{"aggregate-function", required_argument, 0, 'f' },
		{"reduce-function",    required_argument, 0, 'r' },
		{"threads",     optional_argument, 0, 't' },
		{0,0,0,0}
	};

	tomahawk::twk_two_settings settings;
	int32_t x_bins = 1000, y_bins = 1000;
	std::string aggregate_func_name, reduce_func_name;
	int32_t min_cutoff = 5;
	settings.out = "-";
	settings.out_type = 'b';

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:o:O:I:x:y:f:r:c:t:", long_options, &long_index)) != -1){
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
		case 'o':
			settings.out = std::string(optarg);
			break;
		case 'O':
			settings.out_type = std::string(optarg)[0];
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

	if(settings.out_type != 'u' && settings.out_type != 'b'){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Output format must be either 'u' or 'b'..." << std::endl;
		return(1);
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

	std::ofstream* writer; bool use_writer = false;
	if(settings.out.size() == 1 && settings.out == "-"){

	} else {
		use_writer = true;
		writer = new std::ofstream(settings.out, std::ios::out);
		if(!writer->good()){
			std::cerr << tomahawk::utility::timestamp("ERROR") << "Could not open output file \"" << settings.out << "\"!" << std::endl;
			return(1);
		}

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
		else if(reduce_func_name == "count"){ r = &tomahawk::twk_sstats::GetCount; }
		else if(reduce_func_name == "n")    { r = &tomahawk::twk_sstats::GetCount; }
		else if(reduce_func_name == "total"){ r = &tomahawk::twk_sstats::GetTotal; }
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
    if(oreader.BuildIntervals(settings.ivals,oreader.hdr.GetNumberContigs(),oreader.index,oreader.hdr) == false){
        return 1;
    }

    // If we have settings then we do not need to perform two-pass over data as
    // we already know the scene dimensions.
	if(settings.ivals.size()){

	}

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
	uint32_t n_chrom_set = 0;
	for(int i = 0; i < slaves[0].contig_avail.size(); ++i){
		n_chrom_set += slaves[0].contig_avail[i].set;
	}

	// Step 2: Determine boundaries given the contigs that were set.
	//         Calculate the landscape ranges (X and Y dimensions).
	// Approach 2: Dropping regions with no data.
	uint64_t range = 0;
	std::vector<tomahawk::twk1_aggregate_t::offset_tuple> rid_offsets(slaves[0].contig_avail.size());

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

	tomahawk::twk1_aggregate_t agg(x_bins, y_bins);
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
	slaves[0].Overload(agg, min_cutoff);

	if(use_writer){
		if(settings.out_type == 'b'){
			*writer << agg;
		} else {
			slaves[0].PrintMatrix(*writer, min_cutoff);
		}
		writer->flush();
		writer->close();
		delete writer;
	} else {
		if(settings.out_type == 'b') std::cout << agg;
		else slaves[0].PrintMatrix(std::cout, min_cutoff);
	}

	std::cerr << tomahawk::utility::timestamp("LOG") << "Aggregated " << tomahawk::utility::ToPrettyString(n_recs) << " records in " << tomahawk::utility::ToPrettyString(x_bins*y_bins) << " bins." << std::endl;
	std::cerr << tomahawk::utility::timestamp("LOG") << "Finished." << std::endl;

	delete[] slaves;

	return(0);
}
