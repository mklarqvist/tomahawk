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

	tomahawk::twk_two_writer_t writer;

	writer.oindex.SetChroms(oreader.hdr.GetNumberContigs());
	writer.Open(out); writer.mode = 'b';
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

	std::cerr << "testing" << std::endl;
	std::ifstream re(out, std::ios::binary | std::ios::in | std::ios::ate);
	if(re.good() == false){
		std::cerr << "failed reopen" << std::endl;
		return false;
	}
	const uint64_t filesize = re.tellg();
	//re.seekg(0);

	obuf.reset(); obuf2.reset();
	obuf.resize(256000);
	obuf2.resize(256000);

	for(int p = 0; p < offsets.size(); ++p){
		std::cerr << "seeking to: " << offsets[p].foff << std::endl;
		re.seekg(offsets[p].foff);
		if(re.good() == false){
			std::cerr << "failed seek" << std::endl;
			return 1;
		}

		size_t left = offsets[p].nc;
		size_t n_uncomp = 0;
		oreader.zcodec.InitStreamDecompress();
		uint32_t iterations = 0;
		while(left){
			//if(iterations++ == 2) return 1;
			if(re.good() == false){
				std::cerr << "stream is dead" << std::endl;
				return 1;
			}

			// Read a block of data
			std::cerr << "reader=" << re.tellg() << "/" << filesize << " left=" << left << std::endl;
			const size_t n_read = left < 256000 ? left : 256000;
			re.read(obuf.data() + obuf.n_chars_, n_read);
			if(re.good() == false){
				std::cerr << "failed to read: " << n_read << std::endl;
				return 1;
			}
			obuf.n_chars_ += n_read;

			// Decompress
			std::cerr << "read=" << n_read << " and " << obuf.iterator_position_ << "/" << obuf.n_chars_ << std::endl;
			if(oreader.zcodec.StreamDecompress(obuf, obuf2) == false){
				std::cerr << "fail stream decompress" << std::endl;
				return 1;
			}

			// Print head records
			/*
			uint32_t n_rcs = obuf2.size() / tomahawk::twk1_two_t::packed_size;
			tomahawk::twk1_two_t rec;
			obuf2.iterator_position_ = 0;
			for(int i = 0; i < 5; ++i){
				std::cerr << "HEAD\t"<< i << "/" << n_rcs << "@" << obuf2.iterator_position_ << "/" << obuf2.size() << " ";
				obuf2 >> rec;
				rec.PrintLD(std::cerr);
			}
			obuf2.iterator_position_ = (tomahawk::twk1_two_t::packed_size * (n_rcs - 5));
			std::cerr << "moving up to: " << obuf2.iterator_position_ << "/" << obuf2.size() << std::endl;


			std::cerr << "records=" << n_rcs << std::endl;

			for(int i = n_rcs - 5; i < n_rcs; ++i){
				std::cerr << "TAIL\t" << i << "/" << n_rcs << "@" << obuf2.iterator_position_ << "/" << obuf2.size() << " ";
				obuf2 >> rec;
				rec.PrintLD(std::cerr);
			}
*/
			// If input data is not all consumed
			std::cerr << "decomp=" << obuf.iterator_position_ << " -> " << obuf2.size() << std::endl;
			n_uncomp += obuf2.size();
			if(obuf.n_chars_ - obuf.iterator_position_ != 0){
				std::cerr << "moving=" << obuf.n_chars_ - obuf.iterator_position_ << "(" << obuf.n_chars_ << " and " << obuf.iterator_position_ << ")" << std::endl;
				std::memmove(obuf.data(), obuf.data() + obuf.iterator_position_, obuf.n_chars_ - obuf.iterator_position_);
				obuf.n_chars_ = obuf.n_chars_ - obuf.iterator_position_;
			} else obuf.n_chars_ = 0;
			obuf.iterator_position_ = 0;

			//std::cerr << "division=" << (float)obuf2.size() / tomahawk::twk1_two_t::packed_size << " size=" << tomahawk::twk1_two_t::packed_size << std::endl;
			const uint32_t residual = obuf2.size() % tomahawk::twk1_two_t::packed_size;

			//std::cerr << "move=" << obuf2.size() / tomahawk::twk1_two_t::packed_size * tomahawk::twk1_two_t::packed_size << "/" << obuf2.size() << " total of=" << residual << std::endl;
			std::memmove(obuf2.data(), obuf2.data() + (obuf2.size() - residual), residual);
			obuf2.n_chars_ = residual;
			obuf2.iterator_position_ = 0;
			n_uncomp -= residual;

			left -= n_read;
			assert(re.good());
			//obuf.reset();
		}
		assert(n_uncomp == offsets[p].n);
		//std::cerr << "left=" << left << " and " << re.tellg() << "/" << offsets[p].fend << std::endl;
		//std::cerr << "done reading=" << offsets[p].nc << "->" << n_uncomp << "/" << offsets[p].n << std::endl;
		//std::cerr << "remainder=" << oreader.zcodec.outbuf.pos << "/" << oreader.zcodec.outbuf.size << std::endl;
		//std::cerr << "remainder-in=" << oreader.zcodec.inbuf.pos << "/" << oreader.zcodec.inbuf.size << std::endl;
		//std::cerr << "expect=" << offsets[p].n << " and " << offsets[p].nc << " adding " << n_uncomp+oreader.zcodec.outbuf.pos << std::endl;
	}

	// destream

	return 0;
}
