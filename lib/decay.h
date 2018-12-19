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

void decay_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Compute LD decay over distance.\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " decay [options] -i <in.two>\n\n"
	"Options:\n"
	"  -i FILE   input TWO file (required)\n"
	"  -I STRING interval string for target region\n"
	"  -m        output haplotypes in tab-delimited matrix form\n\n";
}

int decay(int argc, char** argv){
	if(argc < 3){
		decay_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{"intervals",   required_argument, 0, 'I' },
		{"numeric",     optional_argument, 0, 'n' },
		{"matrix",      no_argument,       0, 'm' },
		{0,0,0,0}
	};

	tomahawk::twk_two_settings settings;
	std::vector<std::string> intervals;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:I:?", long_options, &long_index)) != -1){
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
		case 'I':
			intervals.push_back(std::string(optarg));
			break;
		}
	}

	if(settings.in.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(intervals.size() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No interval(s) provided..." << std::endl;
		return(1);
	}

	// Print messages
	tomahawk::ProgramMessage();
	std::cerr << tomahawk::utility::timestamp("LOG") << "Calling decay..." << std::endl;

	// New instance of reader.
	tomahawk::two_reader oreader;

	// Open file handle.
	if(oreader.Open(settings.in) == false){
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

	// Build intervals data structures if any are available.
	if(settings.intervals.Build(settings.ivals,oreader.hdr.GetNumberContigs(),oreader.index,oreader.hdr) == false)
		return 1;

	uint64_t n_range = 10e6;
	uint32_t n_bins  = 1000;
	uint32_t n_range_bin = 10e6/1000;
	std::vector<std::pair<double,uint64_t>> decay(n_range_bin+1,{0,0});

	while(oreader.NextRecord()){
		// Same contig only.
		if(oreader.it.rcd->ridA == oreader.it.rcd->ridB){
			// Upper trig only.
			if(oreader.it.rcd->Apos < oreader.it.rcd->Bpos){
				decay[std::min((oreader.it.rcd->Bpos - oreader.it.rcd->Apos) / n_range_bin, n_bins)].first += oreader.it.rcd->R2;
				++decay[std::min((oreader.it.rcd->Bpos - oreader.it.rcd->Apos) / n_range_bin, n_bins)].second;
			}
		}
	}

	std::cout << "From\tTo\tMean\tFrequency\n";
	for(int i = 0; i < decay.size(); ++i){
		std::cout << (i*n_range_bin) << '\t' << ((i+1)*n_range_bin) << '\t' << decay[i].first/std::max(decay[i].second,(uint64_t)1) << '\t' << decay[i].second << '\n';
	}
	std::cout.flush();

	return 0;
}
