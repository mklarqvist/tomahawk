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
#include <getopt.h>

#include "tomahawk/two/TomahawkOutputFilterController.h"
#include "tomahawk/two/TomahawkOutputReader.h"
#include "utility.h"
#include "totempole/TotempoleReader.h"
#include "tomahawk/TomahawkReader.h"

void stats_usage(void){
	programMessage();
	std::cerr <<
	"About:  Calculates basic summary statistics for a TWK/TWO file.\n"
	"        Data does not have to be indexed. However, operations are faster if they\n"
	"        are.\n"
	"Usage:  " << Tomahawk::Constants::PROGRAM_NAME << " stats [options] -i <in.two>\n\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -b INT   number of bins (default: 10)\n";
}

int stats(int argc, char** argv){
	if(argc < 3){
		stats_usage();
		return(1);
	}

	static struct option long_options[] = {
		{"input",		required_argument, 0, 'i' },
		{"bins",		optional_argument, 0, 'b' },
		{"silent",		no_argument, 0, 's' },
		{0,0,0,0}
	};

	// Parameter defaults
	std::string input, output;

	int c = 0;
	int long_index = 0;
	S32 bins = 10;
	while ((c = getopt_long(argc, argv, "i:b:s", long_options, &long_index)) != -1){
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
		case 'b':
			bins = atoi(optarg);
			if(bins <= 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Cannot set bins to <= 0..." << std::endl;
				return(1);
			}
			break;
		case 's':
			SILENT = 1;
			break;
		}
	}

	if(input.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No input file specified..." << std::endl;
		return(1);
	}

	if(!SILENT){
		programMessage();
		std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling stats..." << std::endl;
	}

	// Todo: move out
	std::vector<std::string> inputFile_parts = Tomahawk::Helpers::split(input, '.');
	std::string& end = inputFile_parts[inputFile_parts.size() - 1];
	std::transform(end.begin(), end.end(), end.begin(), ::tolower); // transform chars to lower case

	if(end == Tomahawk::Constants::OUTPUT_SUFFIX){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Unsupported" << std::endl;
		return(1);

	} else if(end == Tomahawk::Constants::OUTPUT_LD_SUFFIX){
		Tomahawk::IO::TomahawkOutputReader reader;
		//reader.setWriteHeader(outputHeader);
		//Tomahawk::TomahawkOutputFilterController& filter = reader.getFilter();
		//filter = Tomahawk::TomahawkOutputFilterController(two_filter); // use copy ctor to transfer data

		//if(!reader.setWriterType(outputType))
		//	return 1;

		if(!reader.Open(input))
			return 1;

		//if(!reader.AddRegions(filter_regions)){
		//	std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Failed to add region!" << std::endl;
		//	return 1;
		//}

		if(!reader.summary(input, bins))
			return 1;

	} else {
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Unrecognised input file format: " << input << std::endl;
		return 1;
	}

	return 0;
}
