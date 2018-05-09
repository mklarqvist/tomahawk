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
#include <vector>

#include "tomahawk/two/TomahawkOutputReader.h"
#include "utility.h"

void concat_usage(void){
	programMessage();
	std::cerr <<
	"About:  Concatenate or combine TWO files. All source files must share the same contig\n"
	"        information in the same order as described in the header. This program is used\n"
	"        primarily to concatenate output TWO files from calc in chunk mode.\n"
	"        The input files does not have to be sorted.\n"
	"Usage:  " << Tomahawk::Constants::PROGRAM_NAME << " concat [options] <in.two>\n\n"
	"Options:\n"
	"  -i FILE  input files (required)\n"
	"  -F LIST  list of files to concatenate (required)\n"
	"  -o FILE  output file (required)\n"
	"  -s       Hide all program messages [null]\n";
}

int concat(int argc, char** argv){
	if(argc < 3){
		concat_usage();
		return(1);
	}

	static struct option long_options[] = {
		{"input",		optional_argument, 0, 'i' },
		{"output",		required_argument, 0, 'o' },
		{"files",		optional_argument, 0, 'F' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	// Parameter defaults
	std::vector<std::string> input;
	std::string output, files;

	int c = 0;
	int long_index = 0;
	while ((c = getopt_long(argc, argv, "i:o:F:s", long_options, &long_index)) != -1){
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
			input.push_back(std::string(optarg));
			break;
		case 'o':
			output = std::string(optarg);
			break;
		case 'F':
			files = std::string(optarg);
			break;
		case 's':
			SILENT = 1;
			break;
		}
	}

	if(input.size() == 0 && files.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No input file specified..." << std::endl;
		std::cerr << input.size() << '\t' << files.size() << std::endl;
		return(1);
	}

	if(output.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No output file specified..." << std::endl;
		return(1);
	}

	if(files.length() != 0 && input.size() != 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Cannot specify both list of input files and manually declare input files..." << std::endl;
		return(1);
	}

	if(!SILENT){
		programMessage();
		std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling concat..." << std::endl;
	}

	Tomahawk::TomahawkOutputReader reader;
	if(input.size() == 0){
		if(!reader.concat(files, output)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "CONCAT") << "Failed to concat files!" << std::endl;
			return 1;
		}

	} else {
		if(!reader.concat(input, output)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "CONCAT") << "Failed to concat files!" << std::endl;
			return 1;
		}

	}

	return 0;
}
