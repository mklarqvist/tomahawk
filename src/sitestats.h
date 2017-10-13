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

#include "utility.h"
#include "totempole/TotempoleReader.h"
#include "tomahawk/TomahawkCalculations.h"

void sitestats_usage(void){
	programMessage();
	std::cerr <<
	"About:  Calculates Tajima's D, mean nucleotide diversity, mean minor allele frequency, mean het freq\n"
	"        from Tomahawk file. If groups are defined (-g) then all statistics will be computed\n"
	"        for each group\n"
	"Usage:  " << Tomahawk::Constants::PROGRAM_NAME << " tajima [options] -i <in.twk>\n\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -g FILE  input groups file\n"
	"  -b INT   bin-size in bases (default: 1000)\n"
	"  -B       do not bin (output by chromosome)\n";
}

int sitestats(int argc, char** argv){
	if(argc < 3){
		tajima_usage();
		return(1);
	}

	static struct option long_options[] = {
		{"input",		required_argument, 0, 'i' },
		{"silent",		no_argument, 0, 's' },
		{"groups",		optional_argument, 0, 'g' },
		{0,0,0,0}
	};

	// Parameter defaults
	std::string input, output, groups;

	int c = 0;
	int long_index = 0;
	S32 bin_size = 1000;
	while ((c = getopt_long(argc, argv, "i:g:s", long_options, &long_index)) != -1){
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

		case 'g':
			groups = std::string(optarg);
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
		std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling sitestats..." << std::endl;
	}

	Tomahawk::TomahawkCalculations tomahawk;
	if(!tomahawk.Open(input)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Failed build!" << std::endl;
		return 1;
	}

	if(groups.size() > 0){
	//const std::string temp = "/media/klarqv01/NVMe/1kgp3/populations/integrated_call_samples_v3.20130502.ALL.panel";
		if(!tomahawk.loadGroups(groups)){
			std::cerr << "could not load groups" << std::endl;
			return 1;
		}
	}

	if(!tomahawk.calculateSiteStats()){
		std::cerr << "failed" << std::endl;
		return 1;
	}

	return 0;
}
