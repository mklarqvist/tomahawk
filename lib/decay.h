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
	"  -w INT    window size in base bairs (default: 10Mb)\n"
	"  -b INT    number of bins each window is separated into\n"
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
		{"window",      optional_argument, 0, 'w' },
		{"bins",        optional_argument, 0, 'b' },
		{"matrix",      no_argument,       0, 'm' },
		{0,0,0,0}
	};

	tomahawk::twk_two_settings settings;

	int64_t n_range = 10e6;
	int32_t n_bins  = 1000;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:I:w:b:?", long_options, &long_index)) != -1){
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
			settings.ivals.push_back(std::string(optarg));
			break;
		case 'w':
			n_range = atoi(optarg); break;
		case 'b':
			n_bins = atoi(optarg); break;
		}
	}

	if(settings.in.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(settings.ivals.size() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No interval(s) provided..." << std::endl;
		return(1);
	}

	// Print messages
	tomahawk::ProgramMessage();
	std::cerr << tomahawk::utility::timestamp("LOG") << "Calling decay..." << std::endl;

	// New instance of reader.
	tomahawk::two_reader oreader;
	//if(oreader.Decay(settings, n_range, n_bins) == false) return 1;
	if(oreader.PositionalDecay(settings) == false) return 1;
	return 0;
}
