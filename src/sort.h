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
#include "algorithm/sort/output_sorter.h"
#include "tomahawk/tomahawk_reader.h"
#include "tomahawk/two/tomahawk_output_reader.h"
#include "utility.h"

void sort_usage(void){
	programMessage();
	std::cerr <<
	"About:  Sort TWO files: provides two basic subroutines. Uncompressed files are generally too big to\n"
	"        be sorted in available memory. Because of this, sorting is split into two \n"
	"        subroutines: parallel-partial sort and a merge-sort. First sort the file without -M\n"
	"        triggered and then run sort on that output with -M flag to perform a k-way merge sort\n"
	"        using the partially block-sorted data.\n"
	"        Note that combining -L and -t incur at least O(L*t) memory!\n"
	"Usage:  " << tomahawk::constants::PROGRAM_NAME << " sort [options] <in.two>\n\n"
	"Options:\n"
	"  -i FILE   input Tomahawk (required)\n"
	"  -o FILE   output file (required)\n"
	"  -L FLOAT  memory limit in MB (default: 1e6)\n"
	"  -t INT    threads (default: " + std::to_string(std::thread::hardware_concurrency()) + ")\n"
	"  -M        merge [null]\n"
	"  -b INT    block size in MB when merging (default: 1e6)\n"
	"  -s        Hide all program messages [null]\n";
}

int sort(int argc, char** argv){
	if(argc < 3){
		sort_usage();
		return(1);
	}

	static struct option long_options[] = {
		{"input",		required_argument, 0, 'i' },
		{"output",		required_argument, 0, 'o' },
		{"memory",		optional_argument, 0, 'L' },
		{"threads",		optional_argument, 0, 't' },
		{"merge",		no_argument, 0, 'M' },
		{"block-size",	optional_argument, 0, 'b' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	// Parameter defaults
	std::string input, output;
	double memory_limit = 1e6;
	double block_size   = 1e6;
	bool merge = false;
	int threads = std::thread::hardware_concurrency();

	int c = 0;
	int long_index = 0;
	while ((c = getopt_long(argc, argv, "i:o:L:t:b:dDMs", long_options, &long_index)) != -1){
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
			input = std::string(optarg);
			break;
		case 'o':
			output = std::string(optarg);
			break;
		case 'L':
			memory_limit = atof(optarg);
			if(memory_limit <= 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter L cannot be negative" << std::endl;
				return(1);
			}
			break;
		case 'b':
			block_size = atoi(optarg);
			if(block_size <= 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter b cannot be negative" << std::endl;
				return(1);
			}
			break;
		case 't':
			threads = atoi(optarg);
			if(threads <= 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter t cannot be <= 0" << std::endl;
				return(1);
			}
			break;

		case 'M':
			merge = true;
			break;
		}
	}

	if(input.length() == 0){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "No input file specified..." << std::endl;
		std::cerr << input.size() << '\t' << input << std::endl;
		return(1);
	}

	if(output.length() == 0){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "No output file specified..." << std::endl;
		std::cerr << output.size() << '\t' << input << std::endl;
		return(1);
	}

	if(!SILENT){
		programMessage();
		std::cerr << tomahawk::helpers::timestamp("LOG") << "Calling sort..." << std::endl;
	}

	tomahawk::algorithm::OutputSorter reader;
	reader.n_threads = threads;

	if(merge == false){
		if(!reader.sort(input, output, memory_limit)){
			std::cerr << tomahawk::helpers::timestamp("ERROR", "SORT") << "Failed to sort file!" << std::endl;
			return 1;
		}
	} else {
		if(!reader.sortMerge(input, output, block_size)){
			std::cerr << tomahawk::helpers::timestamp("ERROR", "SORT") << "Failed merge" << std::endl;
			return 1;
		}
	}

	return 0;
}
