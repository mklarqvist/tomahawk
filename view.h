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

void view_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Convert binary TWO->LD/TWO, subset and slice TWO data\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " view [options] -i <in.two>\n\n"
	"Options:\n"
	"  -i FILE   input TWO file (required)\n"
	"  -h/H      (twk/two) header only / no header\n"
	"  -I STRING filter interval <contig>:pos-pos (TWK/TWO) or linked interval <contig>:pos-pos,<contig>:pos-pos (TWO only)\n\n"
	//"  -J        output JSON object\n\n"
	"  -o FILE    output file (- for stdout; default: -)\n"
	"  -O <b|u>   b: compressed TWK, u: uncompressed VCF\n\n"


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

int view(int argc, char** argv){
	if(argc < 3){
		view_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{"output",      optional_argument, 0, 'o' },
		{"output-type", optional_argument, 0, 'O' },
		{"minP",        optional_argument, 0, 'p' },
		{"maxP",        optional_argument, 0, 'P' },
		{"minR",        optional_argument, 0, 'z' },
		{"maxR",        optional_argument, 0, 'Z' },
		{"minR2",       optional_argument, 0, 'r' },
		{"maxR2",       optional_argument, 0, 'R' },
		{"minDP",       optional_argument, 0, 'b' },
		{"maxDP",       optional_argument, 0, 'B' },
		{"minD",        optional_argument, 0, 'd' },
		{"maxD",        optional_argument, 0, 'D' },
		{"minP1",  optional_argument, 0, '1' },
		{"minP2",  optional_argument, 0, '2' },
		{"minQ1",  optional_argument, 0, '3' },
		{"minQ2",  optional_argument, 0, '4' },
		{"maxP1",  optional_argument, 0, '5' },
		{"maxP2",  optional_argument, 0, '6' },
		{"maxQ1",  optional_argument, 0, '7' },
		{"maxQ2",  optional_argument, 0, '8' },
		{"minMHC", optional_argument, 0, 'a' },
		{"maxMHC", optional_argument, 0, 'A' },

		{"minChi",      optional_argument, 0, 'x' },
		{"maxChi",      optional_argument, 0, 'X' },
		{"minMCV",      optional_argument, 0, 'm' },
		{"maxMCV",      optional_argument, 0, 'M' },

		//{"JSON",        optional_argument, 0, 'J' },
		{"flagInclude", optional_argument, 0, 'f' },
		{"flagExclude", optional_argument, 0, 'F' },
		{"upperTriangular", no_argument, 0, 'u' },
		{"lowerTriangular", no_argument, 0, 'l' },

		{"headerOnly",  no_argument, 0, 'H' },
		{"noHeader",    no_argument, 0, 'h' },

		{"interval",    optional_argument, 0, 'I' },

		{0,0,0,0}
	};

	std::string in, out;
	tomahawk::twk_two_filter filter;
	tomahawk::twk_two_writer_t writer;
	bool write_header = true, header_only = false;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:HhI:o:O:r:R:z:Z:p:P:d:D:b:B:1:2:3:4:5:6:7:8:x:X:a:A:m:M:f:F:ul", long_options, &long_index)) != -1){
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
		case 'O':
			if(std::string(optarg).size() != 1){
				std::cerr << "illegal O" << std::endl;
				return 1;
			}

			writer.mode = optarg[0];
			break;

		case 'p':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetPLow(atof(optarg));
			break;
		case 'P':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetPHigh(atof(optarg));
			break;

		case 'z':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetRLow(atof(optarg));
			break;
		case 'Z':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetRHigh(atof(optarg));
			break;
		case 'r':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetR2Low(atof(optarg));
			break;
		case 'R':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetR2High(atof(optarg));
			break;

		case 'b':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetDprimeLow(atof(optarg));
			break;
		case 'B':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetDprimeHigh(atof(optarg));
			break;

		case 'd':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetDLow(atof(optarg));
			break;
		case 'D':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetDHigh(atof(optarg));
			break;

		case '1':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetHapALow(atof(optarg));
			break;

		case '2':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetHapBLow(atof(optarg));
			break;

		case '3':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetHapCLow(atof(optarg));
			break;

		case '4':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetHapDLow(atof(optarg));
			break;

		case '5':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetHapAHigh(atof(optarg));
			break;

		case '6':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetHapBHigh(atof(optarg));
			break;

		case '7':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetHapCHigh(atof(optarg));
			break;

		case '8':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetHapDHigh(atof(optarg));
			break;

		case 'a':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetMHCLow(atof(optarg));
			break;

		case 'A':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetMHCHigh(atof(optarg));
			break;

		case 'f':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_NUMBER) == false){
				std::cerr << "not a valid number" << std::endl;
				return 1;
			}
			filter.SetFlagInclude(atof(optarg));
			break;

		case 'F':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_NUMBER) == false){
				std::cerr << "not a valid number" << std::endl;
				return 1;
			}
			filter.SetFlagExclude(atof(optarg));
			break;

		case 'x':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetChiSqLow(atof(optarg));
			break;

		case 'X':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetChiSqHigh(atof(optarg));
			break;

		case 'm':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetChiSqModelLow(atof(optarg));
			break;

		case 'M':
			if(std::regex_match(std::string(optarg), tomahawk::TWK_REGEX_FLOATING_EXP) == false){
				std::cerr << "not a valid float" << std::endl;
				return 1;
			}
			filter.SetChiSqModelHigh(atof(optarg));
			break;

		case 'u': filter.SetUpperTrig(); break;
		case 'l': filter.SetLowerTrig(); break;
		case 'h': header_only = true; break;
		case 'H': write_header = false; break;
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

	std::string view_string = "\n##tomahawk_viewVersion=" + std::to_string(VERSION) + "\n";
	view_string += "##tomahawk_viewCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + tomahawk::utility::datetime(); + "\n";
	oreader.hdr.literals_ += view_string;

	writer.oindex.SetChroms(oreader.hdr.GetNumberContigs());
	writer.Open(out);
	if(writer.mode == 'u' && header_only){
		writer.WriteHeader(oreader);
		return 0;
	} else if(writer.mode == 'u' && write_header)
		writer.WriteHeader(oreader);
	else if(writer.mode == 'u' && write_header == false){
		std::cout << "flags\tridA\tposA\tridB\tposB\tHOMHOM\tHOMALT\tALTHOM\tALTALT\tD\tDprime\tR\tR2\tP\tChiSqFisher\tChiSqModel" << std::endl;
	}
	else if(writer.mode == 'b')
		writer.WriteHeader(oreader);

	//tomahawk::twk1_two_iterator it;
	//it.stream = oreader.stream;
	tomahawk::twk1_two_block_t blk2;
	//tomahawk::ZSTDCodec zcodec;
	//std::cerr << "resizing to=" << 1000000000/sizeof(tomahawk::twk1_two_t) << " entries..." << std::endl;
	//blk2.resize(1000000000/sizeof(tomahawk::twk1_two_t));

	//std::FILE* tmpf = std::tmpfile();
	tomahawk::twk_buffer_t obuf;

	//tomahawk::twk_two_filter filter;
	filter.Build();

	//std::cout << "{\"data\":[";
	std::cerr << "before first block" << std::endl;
	while(oreader.NextRecord()){
		if(filter.Filter(oreader.it.rcd)){
			// get first
			break;
		}
	}
	writer.Add(*oreader.it.rcd);
	//oreader.it.rcd->PrintLD(std::cout);

	while(oreader.NextRecord()){
		if(filter.Filter(oreader.it.rcd)){
			//std::cout.put(',');
			//oreader.it.rcd->PrintLD(std::cout);
			//assert(oreader.it.rcd->Apos != 0);
			writer.Add(*oreader.it.rcd);
		}
	}

	//std::cout << "]}\n";
	if(writer.mode == 'b') writer.WriteFinal();
	writer.close();
	return(0);

	uint32_t tot = 0;
	while(oreader.NextBlock()){
		//std::cerr << i++ << " -> " << it.GetBlock().n << std::endl;
		tot += oreader.it.GetBlock().n;
		for(int j = 0; j < oreader.it.blk.n; ++j){
			//it.blk[j].Print(std::cerr);
			if(blk2.n == blk2.m){
				//std::cerr << "resseting" << std::endl;
				blk2.Sort();
				obuf << blk2;
				if(oreader.zcodec.Compress(obuf, oreader.it.oblk.bytes, 1) == false){
					std::cerr << "could not compress" << std::endl;
				}
				oreader.it.oblk.n = obuf.size();
				oreader.it.oblk.nc = oreader.it.oblk.bytes.size();
				std::cerr << oreader.it.oblk.n << "->" << oreader.it.oblk.nc << " " << (float)oreader.it.oblk.n/oreader.it.oblk.nc << "-fold" << std::endl;

				for(int k = 0; k < blk2.n; ++k){
					blk2.rcds[k].PrintLD(std::cout);
				}

				obuf.reset();
				blk2.reset();
			}
			blk2 += oreader.it.blk[j];
			//it.blk[j].Print(std::cout);
		}
		//std::cerr << blk2.n << "/" << blk2.m << std::endl;
	}
	if(blk2.n){
		blk2.Sort();
		obuf << blk2;
		if(oreader.zcodec.Compress(obuf, oreader.it.oblk.bytes, 1) == false){
			std::cerr << "could not compress" << std::endl;
		}
		oreader.it.oblk.n = obuf.size();
		oreader.it.oblk.nc = oreader.it.oblk.bytes.size();
		std::cerr << oreader.it.oblk.n << "->" << oreader.it.oblk.nc << " " << (float)oreader.it.oblk.n/oreader.it.oblk.nc << "-fold" << std::endl;


		for(int k = 0; k < blk2.n; ++k){
			blk2.rcds[k].PrintLD(std::cout);
		}

		obuf.reset();
		blk2.reset();
	}
	std::cerr << "total=" << tot << "/" << oreader.index.GetTotalVariants() << std::endl;


	return 0;
}
