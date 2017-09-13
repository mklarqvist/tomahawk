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
#include "tomahawk/TomahawkReader.h"
#include "tomahawk/TomahawkOutput/TomahawkOutputFilterController.h"
#include "tomahawk/TomahawkOutput/TomahawkOutputReader.h"

void view_usage(void){
	programMessage();
	std::cerr <<
	"About:  Convert binary TWK/TWO to natural text (VCF/LD); subset and slice data\n"
	"        Data does not have to be indexed. However, operations are faster if they\n"
	"        are.\n"
	"Usage:  " << Tomahawk::Constants::PROGRAM_NAME << " view [options] <in.two>\n\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -o FILE  output file (- for stdout)\n"
	"  -h/H     (twk/two) header only / no header [null]\n"
	"  -O char  output type: b for TWO format, n for tab-delimited format (default: b)\n"
	"  -N       output in tab-delimited text format (see -O) [null]\n"
	"  -B       output in binary TWO/TWK format (see -O, default)[null]\n"
	"  -t INT   number of CPU threads (default: maximum available)\n\n"

	// Twk parameters
	"Twk parameters\n"
	"  -G       drop genotypes in output [null]\n\n"

	// Two parameters
	"Two parameters\n"
	"  -p FLOAT   smallest P-value (default: 0)\n"
	"  -P FLOAT   largest P-value (default: 1)\n"
	"  -d FLOAT   smallest D value (default: -1)\n"
	"  -D FLOAT   largest D value (default: 1)\n"
	"  -d FLOAT   smallest D' value (default: -1)\n"
	"  -D FLOAT   largest D' value (default: 1)\n"
	"  -x FLOAT   smallest Chi-squared CV (default: 0)\n"
	"  -X FLOAT   largest Chi-squared CV (default: inf)\n"
	"  -a INT     minimum minor-haplotype count (default: 0)\n"
	"  -A INT     maximum minor-haplotype count (default: inf)\n"
	"  -f INT     include FLAG value\n"
	"  -F INT     exclude FLAG value\n"
	"  -r FLOAT   Pearson's R-squared minimum cut-off value (default: 0.1)\n"
	"  -R FLOAT   Pearson's R-squared maximum cut-off value (default: 1.0)\n"
	"  --p1 FLOAT largest Chi-squared CV (default: inf)\n"
	"  --p2 FLOAT largest Chi-squared CV (default: inf)\n"
	"  --q1 FLOAT largest Chi-squared CV (default: inf)\n"
	"  --q2 FLOAT largest Chi-squared CV (default: inf)\n"
	"  --min<cell> FLOAT largest Chi-squared CV (default: inf)\n"
	"  --max<cell> FLOAT largest Chi-squared CV (default: inf)\n"
	"  -d       Show real-time progress update in cerr [null]\n"
	"  -s       Hide all program messages [null]\n";
}

int view(int argc, char** argv){
	if(argc < 3){
		view_usage();
		return(1);
	}

	static struct option long_options[] = {
		{"input",		required_argument, 0, 'i' },
		{"output",		optional_argument, 0, 'o' },
		{"minP",		optional_argument, 0, 'p' },
		{"maxP",		optional_argument, 0, 'P' },
		{"minR2",		optional_argument, 0, 'r' },
		{"maxR2",		optional_argument, 0, 'R' },
		{"minDprime",	optional_argument, 0, 'd' },
		{"maxDprime",	optional_argument, 0, 'D' },
		{"minAlelles",	optional_argument, 0, 'a' },
		{"maxAlleles",	optional_argument, 0, 'A' },
		{"flagInclude",	optional_argument, 0, 'f' },
		{"flagExclude",	optional_argument, 0, 'F' },
		{"minHWE",	optional_argument, 0, 'w' },
		{"maxHWE",	optional_argument, 0, 'W' },
		{"headerOnly",	no_argument, 0, 'H' },
		{"noHeader",	no_argument, 0, 'h' },
		{"dropGenotypes",	optional_argument, 0, 'G' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	// Parameter defaults
	std::string input, output;
	double minR2 = 0, maxR2 = 1;
	double minP = 0, maxP = 1;
	float minDprime = -1, maxDprime = 1;
	int64_t minAlleles = 0, maxAlleles = std::numeric_limits<int64_t>::max();
	U16 flagInclude = 0, flagExclude = 0;
	bool outputHeader = true;
	int outputType = 1;
	bool dropGenotypes = false;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:o:P:p:a:A:R:r:f:F:d:D:w:W:O:hHGsNB", long_options, &long_index)) != -1){
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
			input = std::string(optarg);
			break;
		case 'o':
			output = std::string(optarg);
			break;
		case 'r':
			minR2 = atof(optarg);
			if(minR2 < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter r cannot be negative" << std::endl;
				return(1);
			}
			if(minR2 > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter r has to be in range 0 < r < 1" << std::endl;
				return(1);
			}
			break;
		case 'R':
			maxR2 = atof(optarg);
			if(maxR2 < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter R cannot be negative" << std::endl;
				return(1);
			}
			if(maxR2 > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter R has to be in range 0 < R < 1" << std::endl;
				return(1);
			}
			break;
		case 'd':
			minDprime = atof(optarg);
			if(minDprime < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter d cannot be negative" << std::endl;
				return(1);
			}
			if(minDprime > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter d has to be in range 0 < d < 1" << std::endl;
				return(1);
			}
			break;
		case 'D':
			maxDprime = atof(optarg);
			if(maxDprime < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter D cannot be negative" << std::endl;
				return(1);
			}
			if(maxDprime > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter D has to be in range 0 < D < 1" << std::endl;
				return(1);
			}
			break;
		case 'p':
			minP = atof(optarg);
			if(minP < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter p cannot be negative" << std::endl;
				return(1);
			}
			if(minP > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter p has to be in range 0 < p < 1" << std::endl;
				return(1);
			}
			break;
		case 'P':
			maxP = atof(optarg);
			if(maxP < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter P cannot be negative" << std::endl;
				return(1);
			}
			if(maxP > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter P has to be in range 0 < P < 1" << std::endl;
				return(1);
			}
			break;
		case 'a':
			minAlleles = atoi(optarg);
			if(minAlleles < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter a cannot be negative" << std::endl;
				return(1);
			}
			break;

		case 'A':
			maxAlleles = atoi(optarg);
			if(maxAlleles < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter A cannot be negative" << std::endl;
				return(1);
			}
			break;

		case 'f':
			flagInclude = atoi(optarg);
			break;
		case 'F':
			flagExclude = atoi(optarg);
			break;
		case 's':
			SILENT = 1;
			--hits;
			break;
		case 'h':
			outputHeader = true;
			--hits;
			break;
		case 'H':
			outputHeader = false;
			--hits;
			break;
		case 'N':
			outputType = 1;
			--hits;
			break;
		case 'B':
			outputType = 0;
			--hits;
			break;
		case 'O':
			outputType = atoi(optarg);
			break;

		case 'G':
			dropGenotypes = true;
			break;
		}
	}

	if(input.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No input file specified..." << std::endl;
		std::cerr << input.size() << '\t' << input << std::endl;
		return(1);
	}

	if(!SILENT){
		programMessage();
		std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling view..." << std::endl;
	}

	// Todo: move out
	std::vector<std::string> inputFile_parts = Tomahawk::Helpers::split(input, '.');
	std::string& end = inputFile_parts[inputFile_parts.size() - 1];
	std::transform(end.begin(), end.end(), end.begin(), ::tolower); // transform chars to lower case

	// Todo: action
	// Parse remainder parameters
	// Assume these parameters are contig or position values for filtering
	std::vector<std::string> filter_regions;
	for(U32 i = 2+hits; i < argc; ++i){
		std::string param(&argv[i][0]);

		std::cerr << param << std::endl;

		if(!Tomahawk::Helpers::parsePositionalStringTWO(param)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Illegal parse of input string: " << param << std::endl;
			return 1;
		}

		filter_regions.push_back(param);
	}

	if(end == Tomahawk::Constants::OUTPUT_SUFFIX){
		Tomahawk::TomahawkReader tomahawk;
		tomahawk.setDropGenotypes(dropGenotypes);
		tomahawk.setShowHeader(outputHeader);
		if(!tomahawk.Open(input)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Failed build!" << std::endl;
			return 1;
		}

		//this->SelectWriterOutputType(Tomahawk::IO::GenericWriterInterace::type::cout);
		tomahawk.outputBlocks();

	} else if(end == Tomahawk::Constants::OUTPUT_LD_SUFFIX){
		Tomahawk::IO::TomahawkOutputReader reader;
		Tomahawk::TomahawkOutputFilterController& filter = reader.getFilter();
		// Todo: move into class
		// Set filter parameters
		if(!filter.setFilterRsquared(minR2, maxR2)) return 1;
		filter.setFilterInclude(flagInclude);
		filter.setFilterExclude(flagExclude);
		if(!filter.setFilterMHF(minAlleles, maxAlleles)) return 1;
		if(!filter.setFilterP(minP, maxP)) return 1;
		if(!filter.setFilterDprime(minDprime, maxDprime)) return 1;
		reader.setWriteHeader(outputHeader);

		if(!reader.setWriterType(outputType))
			return 1;

		if(!reader.Open(input))
			return 1;

		if(!reader.AddRegions(filter_regions)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Failed to add region!" << std::endl;
			return 1;
		}

		if(!reader.view(input))
			return 1;

	} else {
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Unrecognised input file format: " << input << std::endl;
		return 1;
	}

	return 0;
}
