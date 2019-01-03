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

template <class int_t>
	struct pair {
		pair() : t(0), n(0){}

		void operator+=(const int_t v){ t += v; ++n; }

		int_t t;
		uint64_t n;
	};

void stats_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Aggregate TWO data into a rasterized matrix of size [x,y] for\n"
    "        plotting. Data has to be sorted and pre-sliced to the correct\n"
	"        interval-of-interest prior to running.\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " aggregate [options] <in.two>\n\n"

	"Options:\n"
	"  -i   FILE   input TWO file (required)\n"
	"  -x,y INT    number of X/Y-axis bins (default: 1000)\n"
	"  -f   STRING aggregation function: can be one of (r2,r,d,dprime,dp,p,hets,alts,het,alt)(required)\n"
	"  -r   STRING reduction function: can be one of (mean,count,n,min,max,sd)(required)\n"
	"  -I   STRING filter interval <contig>:pos-pos (TWK/TWO) or linked interval <contig>:pos-pos,<contig>:pos-pos\n"
	"  -c   INT    min cut-off value used in reduction function: value < c will be set to 0 (default: 5)\n\n";
}

int stats(int argc, char** argv){
	if(argc < 3){
		stats_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{0,0,0,0}
	};

	tomahawk::twk_two_settings settings;
	int32_t bins = 100;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:", long_options, &long_index)) != -1){
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
		}
	}

	if(settings.in.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	// New instance of reader.
	tomahawk::two_reader oreader;

	// Open file handle.
	if(oreader.Open(settings.in) == false) return 1;


	std::vector< pair<double> > r2(101);
	std::vector< uint64_t > stats(16, 0);
	std::vector< uint64_t > h1(2*oreader.hdr.GetNumberSamples(), 0);
	std::vector< uint64_t > h2(2*oreader.hdr.GetNumberSamples(), 0);
	std::vector< uint64_t > h3(2*oreader.hdr.GetNumberSamples(), 0);
	std::vector< uint64_t > h4(2*oreader.hdr.GetNumberSamples(), 0);

	// Contig-contig matrix.
	std::vector< std::vector<uint64_t> > cmatrix(oreader.hdr.GetNumberContigs(), std::vector<uint64_t>());
	for(int i = 0; i < oreader.hdr.GetNumberContigs(); ++i) cmatrix[i].resize(oreader.hdr.GetNumberContigs());

	while(oreader.NextRecord()){
		r2[uint32_t(oreader.it.rcd->R2 * 100)] += oreader.it.rcd->R2;
		++h1[oreader.it.rcd->cnt[0]];
		++h2[oreader.it.rcd->cnt[1]];
		++h3[oreader.it.rcd->cnt[2]];
		++h4[oreader.it.rcd->cnt[3]];
		++cmatrix[oreader.it.rcd->ridA][oreader.it.rcd->ridB];

		for(int j = 0; j < 16; ++j){
			stats[j] += (oreader.it.rcd->controller & (1 << j)) != 0;
		}
	}

	for(int i = 0; i < r2.size(); ++i){
		std::cout << i << "\t" << r2[i].t << "\t" << r2[i].n << '\n';
	}
	for(int i = 0; i < stats.size(); ++i){
		std::cout << i << "\t" << stats[i] << '\n';
	}
	for(int i = 0; i < h1.size(); ++i){
		std::cout << i << "\t" << h1[i] << "\t" << h2[i] << "\t" << h3[i] << "\t" << h4[i] << '\n';
	}

	std::cout << "contig";
	for(int i = 0; i < oreader.hdr.GetNumberContigs(); ++i){
	    std::cout << '\t' << oreader.hdr.GetContig(i)->name;
	}
	 std::cout.put('\n');

	for(int i = 0; i < oreader.hdr.GetNumberContigs(); ++i){
	    std::cout << oreader.hdr.GetContig(i)->name;
	    for(int j = 0; j < oreader.hdr.GetNumberContigs(); ++j){
	        std::cout << '\t' << cmatrix[i][j];
	    }
	    std::cout.put('\n');
	}

	std::cout.flush();

	return 0;
}
