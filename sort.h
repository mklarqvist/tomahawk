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

void sort_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Convert binary TWO->LD/TWO, subset and slice TWO data\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " view [options] -i <in.two>\n\n"
	"Options:\n"
	"  -i FILE   input TWO file (required)\n"
	"  -h/H      (twk/two) header only / no header\n"
	"  -I STRING filter interval <contig>:pos-pos (TWK/TWO) or linked interval <contig>:pos-pos,<contig>:pos-pos\n\n"
	//"  -J        output JSON object\n\n"
	"  -o FILE    output file (- for stdout; default: -)\n"
	"  -O <b|u>   b: compressed TWO, u: uncompressed LD\n\n"


	// Filter parameters
	"Filter parameters:\n"
	"  -r,-R --minR2,--maxR2   FLOAT   Pearson's R-squared min/max cut-off value\n"
	"  -z,-Z --minR,--maxR     FLOAT   Pearson's R min/max cut-off value\n"
	"  -p,-P --minP,--maxP     FLOAT   Min/max P-value (default: [0,1])\n"
	"  -d,-D --minD,--maxD     FLOAT   Min/max D value (default: [-1,1])\n"
	"  -b,-B --minDP,--maxDP   FLOAT   Min/max D' value (default: [0,1])\n"
	"  -1,-5 --minP1,--maxP1   FLOAT   Min/max REF_REF count (default: [0,inf])\n"
	"  -2,-6 --minP2,--maxP2   FLOAT   Min/max REF_ALT count (default: [0,inf])\n"
	"  -3,-7 --minQ1,--maxQ1   FLOAT   Min/max ALT_REF count (default: [0,inf])\n"
	"  -4,-8 --minQ2,--maxQ1   FLOAT   Min/max ALT_ALT count (default: [0,inf])\n"
	"  -a,-A --minMHC,--maxMHC FLOAT   Min/max number of non-major haplotype count (default: [0,inf])\n"
	"  -x,-X --minChi,--maxChi FLOAT   Min/max Chi-squared CV of contingency table (default: [0,inf])\n"
	"  -m,-M --minMCV,--maxMCV FLOAT   Min/max Chi-squared CV of unphased model (default: [0,inf])\n"
	"  -f  INT  include FLAG value\n"
	"  -F  INT  exclude FLAG value\n"
	"  -u       output only the upper triangular values\n"
	"  -l       output only the lower triangular values\n";
}

int sort(int argc, char** argv){
	if(argc < 3){
		sort_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{"output",      optional_argument, 0, 'o' },
		{"memory-usage", optional_argument, 0, 'm' },


		{0,0,0,0}
	};

	std::string in, out;
	int32_t memory_limit = 0;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:o:m:?", long_options, &long_index)) != -1){
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
			in = std::string(optarg);
			break;
		case 'o':
			out = std::string(optarg);
			break;
		case 'm':
			if(std::string(optarg).size() != 1){
				std::cerr << "illegal O" << std::endl;
				return 1;
			}

			memory_limit = atoi(optarg);
			break;


		}
	}

	if(in.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	tomahawk::two_reader oreader;

	if(oreader.Open(in) == false){
		std::cerr << "failed to open" << std::endl;
		return false;
	}

	std::string view_string = "\n##tomahawk_sortVersion=" + std::to_string(VERSION) + "\n";
	view_string += "##tomahawk_sortCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + tomahawk::utility::datetime(); + "\n";

	// generate temp name

	/*std::string BasePath(const std::string& input);
	std::string BaseName(const std::string& input);
	std::string ExtensionName(const std::string& input);
	*/

	std::string extension = tomahawk::twk_two_writer_t::GetExtension(out);
	if(extension != "two"){
		out += ".two";
	}

	// generate temp names
	std::string suffix = tomahawk::twk_two_writer_t::RandomSuffix();
	std::string base_path = tomahawk::twk_two_writer_t::GetBasePath(out);
	std::string base_name = tomahawk::twk_two_writer_t::GetBaseName(out);
	std::string temp_out = (base_path.size() ? base_path + "/" : "") + base_name + "_" + suffix + ".two";

	tomahawk::twk_two_writer_t writer;
	writer.oindex.SetChroms(oreader.hdr.GetNumberContigs());
	std::cerr << tomahawk::utility::timestamp("LOG","WRITER") << "Opening " << (base_name + "_" + suffix + ".two") << std::endl;
	if(writer.Open(temp_out) == false){
		std::cerr << "failed to open=" << temp_out << std::endl;
		return 1;
	}
	writer.mode = 'b'; // force binary mode
	writer.WriteHeader(oreader);

	tomahawk::twk1_two_block_t blk2;
	std::cerr << "resizing to=" << 1000000000/sizeof(tomahawk::twk1_two_t) << " entries..." << std::endl;
	blk2.resize(1000000000/sizeof(tomahawk::twk1_two_t));

	//std::FILE* tmpf = std::tmpfile();
	tomahawk::twk_buffer_t obuf, obuf2;

	struct sort_helper {
		sort_helper() : foff(0), fend(0), n(0), nc(0){}

		uint64_t foff, fend, n, nc;
	};

	std::vector<sort_helper> offsets;

	uint32_t tot = 0;
	while(oreader.NextBlock()){
		//std::cerr << i++ << " -> " << it.GetBlock().n << std::endl;
		tot += oreader.it.GetBlock().n;
		for(int j = 0; j < oreader.it.blk.n; ++j){
			//it.blk[j].Print(std::cerr);
			if(blk2.n == blk2.m){
				//std::cerr << "resseting" << std::endl;
				blk2.Sort();
				//obuf << blk2; // poor memory management:

				sort_helper rec;
				rec.foff = writer.stream.tellp();
				oreader.zcodec.InitStreamCompress(1);
				// Compress chunks of 10k records
				int k = 0;
				for(; k + 10000 < blk2.n; k += 10000){
					for(int l = k; l < k + 10000; ++l) obuf << blk2.rcds[l];
					rec.nc += oreader.zcodec.StreamCompress(obuf, obuf2, writer.stream, tomahawk::twk1_two_t::packed_size * 5000);
					obuf.reset(); obuf2.reset();

				}

				// Compress residual records
				std::cerr << "remainder=" << k << "/" << blk2.n << std::endl;

				for(int l = k; l < blk2.n; ++l) obuf << blk2.rcds[l];
				rec.nc += oreader.zcodec.StreamCompress(obuf, obuf2, writer.stream, tomahawk::twk1_two_t::packed_size * 5000);
				obuf.reset(); obuf2.reset();

				oreader.zcodec.StopStreamCompress();
				std::cerr << "after stop=" << oreader.zcodec.outbuf.pos << "/" << oreader.zcodec.outbuf.size << std::endl;
				writer.stream.write((const char*)oreader.zcodec.outbuf.dst, oreader.zcodec.outbuf.pos);
				rec.nc += oreader.zcodec.outbuf.pos;
				rec.n = blk2.n * tomahawk::twk1_two_t::packed_size;
				rec.fend = writer.stream.tellp();
				offsets.push_back(rec);

				blk2.reset();

			}
			blk2 += oreader.it.blk[j];
			//it.blk[j].Print(std::cout);
		}
		//std::cerr << blk2.n << "/" << blk2.m << std::endl;
	}
	if(blk2.n){
		//std::cerr << "resseting" << std::endl;
		blk2.Sort();
		//obuf << blk2; // poor memory management:

		sort_helper rec;
		rec.foff = writer.stream.tellp();
		oreader.zcodec.InitStreamCompress(1);
		// Compress chunks of 10k records
		int k = 0;
		for(; k + 10000 < blk2.n; k += 10000){
			for(int l = k; l < k + 10000; ++l) obuf << blk2.rcds[l];
			rec.nc += oreader.zcodec.StreamCompress(obuf, obuf2, writer.stream, tomahawk::twk1_two_t::packed_size * 5000);
			obuf.reset(); obuf2.reset();

		}

		// Compress residual records
		std::cerr << "remainder=" << k << "/" << blk2.n << std::endl;

		for(int l = k; l < blk2.n; ++l) obuf << blk2.rcds[l];
		rec.nc += oreader.zcodec.StreamCompress(obuf, obuf2, writer.stream, tomahawk::twk1_two_t::packed_size * 5000);
		obuf.reset(); obuf2.reset();

		oreader.zcodec.StopStreamCompress();
		std::cerr << "after stop=" << oreader.zcodec.outbuf.pos << "/" << oreader.zcodec.outbuf.size << std::endl;
		writer.stream.write((const char*)oreader.zcodec.outbuf.dst, oreader.zcodec.outbuf.pos);
		rec.nc += oreader.zcodec.outbuf.pos;
		rec.n = blk2.n * tomahawk::twk1_two_t::packed_size;
		rec.fend = writer.stream.tellp();
		offsets.push_back(rec);

		blk2.reset();
	}
	writer.flush();
	writer.close();
	std::cerr << "total=" << tot << "/" << oreader.index.GetTotalVariants() << std::endl;

	for(int i = 0; i < offsets.size(); ++i){
		std::cerr << i << ": " << offsets[i].foff << "->" << offsets[i].fend << " " << offsets[i].n << " " << offsets[i].nc << std::endl;
	}

	struct twk_two_stream_iterator {

		bool Open(const std::string file,
				const uint64_t foff,
				const uint64_t fend,
				const uint32_t n_uncompressed,
				const uint32_t n_compressed)
		{
			stream = std::ifstream(file, std::ios::binary);
			if(stream.good() == false){
				std::cerr << "could not open " << file << std::endl;
				return false;
			}

			stream.seekg(foff);
			if(stream.good() == false){
				std::cerr << "seek failed" << std::endl;
			}

			n = 0; it = 0; n_unc_cum = 0;
			b_left = fend - foff;
			it_tot = 0;
			//std::cerr << "data left=" << b_left << ": " << foff << "-" << fend << std::endl;
			n_unc = n_uncompressed; n_cmp = n_compressed;
			n_tot = n_unc / tomahawk::twk1_two_t::packed_size;
			off_start = foff; off_end = fend;
			out.reset(); out2.reset();

			zcodec.InitStreamDecompress();

			return true;
		}

		bool Next(tomahawk::twk1_two_t& rec){
			if(it == n){
				if(NextBlock() == false)
					return false;
			}

			out2 >> rec;
			++it; ++it_tot;
			return true;
		}

		bool NextBlock(uint32_t b_read = 256000){
			if(b_left == 0)
				return false;

			const uint32_t residual = out2.size() % tomahawk::twk1_two_t::packed_size;

			if(residual){
				//std::cerr << "has residual=" << residual << std::endl;
				std::memmove(out2.data(), out2.data() + (out2.size() - residual), residual);
			}

			//b_read = b_read > b_left ? b_left : b_read;
			const size_t n_read = b_left > b_read ? b_read : b_left;
			n = 0; it = 0;

			n_unc_cum += out2.size();
			if(out.n_chars_ - out.iterator_position_ != 0){
				//std::cerr << "moving=" << out.n_chars_ - out.iterator_position_ << "(" << out.n_chars_ << " and " << out.iterator_position_ << ")" << std::endl;
				std::memmove(out.data(), out.data() + out.iterator_position_, out.n_chars_ - out.iterator_position_);
				out.n_chars_ = out.n_chars_ - out.iterator_position_;
			} else out.n_chars_ = 0;

			out.iterator_position_ = 0;
			out2.n_chars_ = residual;
			out2.iterator_position_ = 0;
			n_unc_cum -= residual;

			if(out.capacity() == 0 || b_read > out.capacity()){
				out.resize(b_read);
				out2.resize(b_read);
			}

			if(stream.good() == false){
				std::cerr << "stream is dead" << std::endl;
				return false;
			}

			// Read a block of data
			//std::cerr << "reader=" << stream.tellg() << " left=" << b_left << std::endl;
			stream.read(out.data() + out.n_chars_, n_read);
			if(stream.good() == false){
				std::cerr << "failed to read: " << n_read << std::endl;
				return false;
			}
			out.n_chars_ += n_read;

			// Decompress
			//std::cerr << "read=" << n_read << " and " << out.iterator_position_ << "/" << out.n_chars_ << std::endl;
			if(zcodec.StreamDecompress(out, out2) == false){
				std::cerr << "fail stream decompress" << std::endl;
				return false;
			}

			// If input data is not all consumed
			//std::cerr << "decomp=" << out.iterator_position_ << " -> " << out2.size() << std::endl;

			n = out2.size() / tomahawk::twk1_two_t::packed_size;
			//std::cerr << "n=" << n << " it=" << it << std::endl;

			b_left -= n_read;
			return true;
		}

		uint32_t n, it, n_tot, it_tot; // n_records, it_pos
		uint32_t b_left, n_unc, n_cmp, n_unc_cum;
		uint64_t off_start, off_end;
		std::ifstream stream;
		tomahawk::twk_buffer_t out, out2;
		tomahawk::ZSTDCodec zcodec;
	};

	struct two_queue_entry {
	public:
		two_queue_entry(const tomahawk::twk1_two_t& data, uint32_t streamID) :
			qid(streamID), rec(data)
	    {}

	    inline bool operator<(const two_queue_entry& other) const {
	        return(!(rec < other.rec));
	    }

	public:
	    uint32_t qid;
	    tomahawk::twk1_two_t rec;
	};

	obuf.reset(); obuf2.reset();
	obuf.resize(256000);
	obuf2.resize(256000);

	//uint32_t k = 0;
	std::priority_queue<two_queue_entry> queue;

	twk_two_stream_iterator* its = new twk_two_stream_iterator[offsets.size()];
	tomahawk::twk1_two_t rec;
	uint64_t n_rec_total = 0;
	for(int p = 0; p < offsets.size(); ++p){
		// open iterators
		if(its[p].Open(temp_out, offsets[p].foff, offsets[p].fend, offsets[p].n, offsets[p].nc) == false){
			std::cerr << "failed open" << std::endl;
			return false;
		}

		if(its[p].Next(rec) == false){
			std::cerr << "failed to get next" << std::endl;
			return false;
		}
		std::cerr << "queue-" << p << std::endl;
		//rec.PrintLD(std::cerr);
		queue.push(two_queue_entry(rec, p));

		n_rec_total += offsets[p].n / tomahawk::twk1_two_t::packed_size;
		std::cerr << "total=" << its[p].n_tot << std::endl;
	}

	if(queue.empty()){
		std::cerr << tomahawk::utility::timestamp("ERROR","SORT") << "No data in queue..." << std::endl;
		return false;
	}

	//std::cerr << "top" << std::endl;
	//queue.top().rec.PrintLD(std::cerr);

	tomahawk::twk_two_writer_t owriter;
	owriter.oindex.SetChroms(oreader.hdr.GetNumberContigs());
	std::cerr << "opening=" << out << std::endl;
	if(owriter.Open(out) == false){
		std::cerr << "failed top open" << std::endl;
		return false;
	}
	owriter.mode = 'b';
	owriter.oindex.state = TWK_IDX_SORTED;
	// Write header
	std::string sort_string = "\n##tomahawk_sortVersion=" + std::to_string(VERSION) + "\n";
	sort_string += "##tomahawk_sortCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + tomahawk::utility::datetime(); + "\n";
	oreader.hdr.literals_ += sort_string;
	if(owriter.WriteHeader(oreader) == false){
		std::cerr << "failed to write header" << std::endl;
		return false;
	}


	// Reference
	uint32_t ridA = queue.top().rec.ridA;

	tomahawk::Timer timer; timer.Start();
	uint64_t n_entries_out = 0;
	while(queue.empty() == false){
		// peek at top entry in queue
		const uint32_t id = queue.top().qid;

		if(queue.top().rec.ridA != ridA){
			//std::cerr << "change in ridA=" << queue.top().rec.ridA << "!=" << ridA << " -> " << owriter.oblock.n << std::endl;
			if(owriter.WriteBlock() == false){
				std::cerr << "failed to flush" << std::endl;
				return false;
			}
		}

		//queue.top().rec.PrintLD(std::cout);
		owriter.Add(queue.top().rec);
		ridA = queue.top().rec.ridA;


		++n_entries_out;
		if(n_entries_out % 1000000 == 0)
			std::cerr << tomahawk::utility::timestamp("LOG") << tomahawk::utility::ToPrettyString(n_entries_out) << "/" << tomahawk::utility::ToPrettyString(n_rec_total) << std::endl;

		// remove this record from the queue
		queue.pop();

		while(its[id].Next(rec)){
			if(!(rec < queue.top().rec)){
				queue.push( two_queue_entry(rec, id) );
				break;
			}

			if(rec.ridA != ridA){
				//std::cerr << "change in ridA=" << rec.ridA << "!=" << ridA << " -> " << owriter.oblock.n << std::endl;

				if(owriter.WriteBlock() == false){
					std::cerr << "failed to flush" << std::endl;
					return false;
				}
			}
			//rec.PrintLD(std::cout);
			owriter.Add(rec);
			ridA = rec.ridA;

			++n_entries_out;
			if(n_entries_out % 1000000 == 0)
				std::cerr << tomahawk::utility::timestamp("LOG") << tomahawk::utility::ToPrettyString(n_entries_out) << "/" << tomahawk::utility::ToPrettyString(n_rec_total) << std::endl;
		}
	}

	owriter.flush();
	owriter.WriteFinal();
	owriter.close();
	std::cerr << "merge time=" << timer.ElapsedString() << std::endl;

	std::cerr << "deleting intermediary" << std::endl;

	if( remove( temp_out.c_str() ) != 0 )
	    std::cerr << "Error deleting file=" << out << std::endl;
	  else {
	    std::cerr << "File successfully deleted" << out << std::endl;
	    return 1;
	  }

	delete[] its;
	return 0;
}
