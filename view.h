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
	"About:  Convert binary TWK->VCF or TWO->LD, subset and slice TWK/TWO data\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " view [options] -i <in.two>\n\n"
	"Options:\n"
	"  -i FILE   input Tomahawk (required)\n"
	"  -h/H      (twk/two) header only / no header\n"
	"  -N        output in tab-delimited text format (see -O)\n"
	"  -B        output in binary TWO/TWK format (see -O, default)\n"
	"  -I STRING filter interval <contig>:pos-pos (TWK/TWO) or linked interval <contig>:pos-pos,<contig>:pos-pos (TWO only)\n"
	//"  -J        output JSON object\n"
	"  -s        Hide all program messages\n\n"

	// Twk parameters
	"TWK parameters\n"
	"  -G       drop genotypes in output\n\n"

	// Two parameters
	"TWO parameters\n"
	"  -o FILE   output file (- for stdout; default: -)\n"
	"  -O char   output type: b for TWO format, u for tab-delimited LD format\n"
	"  -r, --minR2  FLOAT   Pearson's R-squared minimum cut-off value\n"
	"  -R, --maxR2  FLOAT   Pearson's R-squared maximum cut-off value\n"
	"  -z, --minR   FLOAT   Pearson's R minimum cut-off value\n"
	"  -Z, --maxR   FLOAT   Pearson's R maximum cut-off value\n"
	"  -p, --minP   FLOAT   smallest P-value (default: 0)\n"
	"  -P, --maxP   FLOAT   largest P-value (default: 1)\n"
	"  -d, --minD   FLOAT   smallest D value (default: -1)\n"
	"  -D, --maxD   FLOAT   largest D value (default: 1)\n"
	"  -b, --minDP  FLOAT   smallest D' value (default: 0)\n"
	"  -B, --maxDP  FLOAT   largest D' value (default: 1)\n"
	"  -x, --minChi FLOAT   smallest Chi-squared CV of contingency table (default: 0)\n"
	"  -X, --maxChi FLOAT   largest Chi-squared CV of contingency table (default: inf)\n"
	"  -a, --minMHC FLOAT   minimum minor-haplotype count (default: 0)\n"
	"  -A, --maxMHC FLOAT   maximum minor-haplotype count (default: inf)\n"
	"  -m, --minMCV FLOAT   smallest Chi-squared CV of unphased model (default: 0)\n"
	"  -M, --maxMCV FLOAT   largest Chi-squared CV of unphased model (default: inf)\n"
	"  -1, --maxP1  FLOAT   maximum REF_REF count\n"
	"  -2, --maxP2  FLOAT   maximum REF_ALT count\n"
	"  -3, --maxQ1  FLOAT   maximum ALT_REF count\n"
	"  -4, --maxQ2  FLOAT   maximum ALT_ALT count\n"
	"  -f           INT     include FLAG value\n"
	"  -F           INT     exclude FLAG value\n"
	"  -u                   output only the upper triangular values\n"
	"  -l                   output only the lower triangular values\n";
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
		{"minChi",      optional_argument, 0, 'x' },
		{"maxChi",      optional_argument, 0, 'X' },
		{"minAlelles",  optional_argument, 0, 'a' },
		{"maxAlleles",  optional_argument, 0, 'A' },
		{"minMCV",      optional_argument, 0, 'm' },
		{"maxMCV",      optional_argument, 0, 'M' },
		{"JSON",        optional_argument, 0, 'J' },
		{"flagInclude", optional_argument, 0, 'f' },
		{"flagExclude", optional_argument, 0, 'F' },
		{"upperTriangular", no_argument, 0, 'u' },
		{"lowerTriangular", no_argument, 0, 'l' },
		{"headerOnly",  no_argument, 0, 'H' },
		{"noHeader",    no_argument, 0, 'h' },
		{"dropGenotypes",optional_argument, 0, 'G' },
		{"interval",    optional_argument, 0, 'I' },
		{"silent",      no_argument, 0,  's' },
		{"maxP1",  optional_argument, 0, '1' },
		{"maxP2",  optional_argument, 0, '2' },
		{"maxQ1",  optional_argument, 0, '3' },
		{"maxQ2",  optional_argument, 0, '4' },

		{0,0,0,0}
	};

	std::string in, out;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:o:", long_options, &long_index)) != -1){
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

	tomahawk::twk1_two_iterator it;
	it.stream = oreader.stream;
	tomahawk::twk1_two_block_t blk2;
	std::cerr << "resizing to=" << 1000000000/sizeof(tomahawk::twk1_two_t) << " entries..." << std::endl;
	blk2.resize(1000000000/sizeof(tomahawk::twk1_two_t));

	//std::FILE* tmpf = std::tmpfile();

	uint32_t tot = 0;
	while(it.NextBlock()){
		//std::cerr << i++ << " -> " << it.GetBlock().n << std::endl;
		tot += it.GetBlock().n;
		for(int j = 0; j < it.blk.n; ++j){
			//it.blk[j].Print(std::cerr);
			if(blk2.n == blk2.m){
				//std::cerr << "resseting" << std::endl;
				blk2.Sort();
				for(int k = 0; k < blk2.n; ++k){
					blk2.rcds[k].Print(std::cout);
				}
				blk2.reset();
			}
			blk2 += it.blk[j];
			//it.blk[j].Print(std::cout);
		}
		//std::cerr << blk2.n << "/" << blk2.m << std::endl;
	}
	if(blk2.n){
		blk2.Sort();
		for(int k = 0; k < blk2.n; ++k){
			blk2.rcds[k].Print(std::cout);
		}
		blk2.reset();
	}
	std::cerr << "total=" << tot << "/" << oreader.index.GetTotalVariants() << std::endl;


	return 0;
}
