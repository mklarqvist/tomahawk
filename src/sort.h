/*
Copyright (C) 2016-2017 Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk21@sanger.ac.uk>

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
#include "utility.h"
#include "tomahawk/TomahawkReader.h"
#include "totempole/TotempoleReader.h"
#include "tomahawk/TomahawkOutput/TomahawkOutputReader.h"
#include "algorithm/sort/TomahawkOutputSort.h"

void sort_usage(void){
	programMessage();
	std::cerr <<
	"About:  Sort TWO files: provides two basic subroutines. If the file is too big to\n"
	"        be sorted in available memory, use the -L option to split the file into\n"
	"        sorted chunks no larger than -L MB in size. Then rerun sort with the -M option\n"
	"        to perform a k-way merge sort using the partially block-sorted data.\n"
	"Usage:  " << Tomahawk::Constants::PROGRAM_NAME << " sort [options] <in.two>\n\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -o FILE  output file (required)\n"
	"  -L INT   memory limit in MB (default: 100)\n"
	"  -M       merge [null]\n"
	"  -s       Hide all program messages [null]\n";
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
		{"merge",		no_argument, 0, 'M' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	// Parameter defaults
	std::string input, output;
	S32 memory_limit = 100e6;
	bool merge = false;

	int c = 0;
	int long_index = 0;
	while ((c = getopt_long(argc, argv, "i:o:L:Ms", long_options, &long_index)) != -1){
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
			memory_limit = atoi(optarg) * 1e6;
			if(memory_limit < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter L cannot be negative" << std::endl;
				return(1);
			}
			break;

		case 'M':
			merge = true;
			break;
		}
	}

	if(input.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No input file specified..." << std::endl;
		std::cerr << input.size() << '\t' << input << std::endl;
		return(1);
	}

	if(output.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No output file specified..." << std::endl;
		std::cerr << output.size() << '\t' << input << std::endl;
		return(1);
	}

	if(!SILENT){
		programMessage();
		std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling sort..." << std::endl;
	}

	Tomahawk::Algorithm::Output::TomahawkOutputSorter reader;
	if(!merge){
		if(!reader.sort(input, output, memory_limit)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "SORT") << "Failed to sort file!" << std::endl;
			return 1;
		}
	} else {
		if(!reader.sortMerge(input, output, 10e6)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "SORT") << "failed merge" << std::endl;
			return 1;
		}
	}

	return 0;
}
