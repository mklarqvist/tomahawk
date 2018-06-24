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
#ifndef AGGREGATE_H_
#define AGGREGATE_H_


#include "tomahawk/two/tomahawk_output_reader.h"
#include "utility.h"

void aggregate_usage(void){
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
	"  -i FILE       input Tomahawk output file (TWO; required)\n"
	"  -x INT        X-axis scene dimension in pixels (required)\n"
	"  -y INT        Y-axis scene dimension in pixels (required)\n"
	"  -X INT-INT    X-axis data range (from, to)\n"
	"  -Y INT-INT    Y-axis data range (from, to)\n"
	"  -d STRING     Data field to aggregate: May be one of {R, R2, D, DPrime, p1, p2, q1, q2, P}(required)\n"
	"  -f STRING     Aggregation function: May be one of {count, min, max, mean, sd}(required)\n"
	"  -s            Hide all program messages [null]\n\n";
}

int aggregate(int argc, char** argv){
	if(argc < 3){
		aggregate_usage();
		return(1);
	}

	static struct option long_options[] = {
		{"input",  required_argument, 0, 'i' },
		{"silent", no_argument,       0, 's' },
		{0,0,0,0}
	};

	// Parameter defaults
	std::string input;

	int c = 0;
	int long_index = 0;
	while ((c = getopt_long(argc, argv, "i:x:y:X:Y:d:f:s", long_options, &long_index)) != -1){
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
		}
	}

	if(input.length() == 0){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "No input file specified..." << std::endl;
		return(1);
	}

	if(!SILENT){
		programMessage();
		std::cerr << tomahawk::helpers::timestamp("LOG") << "Calling stats..." << std::endl;
	}

	tomahawk::TomahawkOutputReader reader;
	if(!reader.open(input))  return 1;
	if(!reader.statistics()) return 1;


	return 0;
}



#endif /* AGGREGATE_H_ */
