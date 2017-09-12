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
#include "tomahawk/TomahawkOutput/TomahawkOutputReader.h"

void concat_usage(void){
	programMessage();
	std::cerr <<
	"Usage: " << Tomahawk::Constants::PROGRAM_NAME << " sort [options] <in.two>\n"
	"\n"
	"Options:\n"
	"  -i FILE  input files (required)\n"
	"  -F LIST  list of files to concatenate (required)\n"
	"  -o FILE  output file (required)\n"
	"  -t INT   number of CPU threads (default: maximum available)\n";
}

int concat(int argc, char** argv){
	if(argc < 3){
		sort_usage();
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
	std::string input, output, files;

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
			input = std::string(optarg);
			break;
		case 'o':
			output = std::string(optarg);
			break;
		case 'F':
			files = std::string(optarg);
			break;
		}
	}

	if(input.length() == 0 && files.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No input file specified..." << std::endl;
		std::cerr << input.size() << '\t' << input << std::endl;
		return(1);
	}

	if(output.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No output file specified..." << std::endl;
		std::cerr << output.size() << '\t' << input << std::endl;
		return(1);
	}

	if(files.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No list of files specified..." << std::endl;
		std::cerr << files.size() << '\t' << input << std::endl;
		return(1);
	}

	if(!SILENT){
		programMessage();
		std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling concat..." << std::endl;
	}

	Tomahawk::IO::TomahawkOutputReader reader;
	if(input.size() == 0) input = files;
	if(!reader.concat(files, output)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "CONCAT") << "Failed to concat files!" << std::endl;
		return 1;
	}

	return 0;
}
