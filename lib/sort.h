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
	"About:  Sort TWO files\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " sort [options] -i <in.two>\n\n"
	"Options:\n"
	"  -i FILE   input TWO file (required)\n"
	"  -o FILE   output file (- for stdout; default: -)\n"
	"  -m FLOAT  maximum memory usage per thread in GB (default: 0.5)\n"
	"  -c INT    compression level 1-20 (default: 1)\n"
	"  -t INT    number of threads (default: maximum available)\n\n";
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
		{"compression-level", optional_argument, 0, 'c' },
		{"threads", optional_argument, 0, 't' },
		{0,0,0,0}
	};

	tomahawk::two_sorter_settings settings;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:o:m:c:t:?", long_options, &long_index)) != -1){
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
		case 'm':
			settings.memory_limit = atof(optarg);
			break;
		case 'c':
			settings.c_level = atoi(optarg);
			break;
		case 't':
			settings.n_threads = atoi(optarg);
			break;
		}
	}

	if(settings.in.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(settings.memory_limit <= 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot set memory limit <= 0..." << std::endl;
		return(1);
	}

	if(settings.n_threads <= 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot set number of threads <= 0..." << std::endl;
		return(1);
	}

	if(settings.c_level <= 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot set the compression level <= 0..." << std::endl;
		return(1);
	}

	// Print messages
	tomahawk::ProgramMessage();
	std::cerr << tomahawk::utility::timestamp("LOG") << "Calling sort..." << std::endl;

	tomahawk::two_reader sorter;
	if(sorter.Sort(settings) == false) return 1;
	return 0;
}
