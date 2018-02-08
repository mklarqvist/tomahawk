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

#include "tomahawk/two/output_filter.h"
#include "tomahawk/two/TomahawkOutputReader.h"
#include "utility.h"
#include "totempole/TotempoleReader.h"
#include "tomahawk/TomahawkReader.h"

void view_usage(void){
	programMessage();
	std::cerr <<
	"About:  Convert binary TWK->VCF or TWO->LD; subset and slice TWK/TWO data\n"
	"        Data does not have to be indexed. However, operations are faster if they\n"
	"        are.\n"
	"Usage:  " << Tomahawk::Constants::PROGRAM_NAME << " view [options] -i <in.two>\n\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -o FILE  output file (- for stdout)\n"
	"  -h/H     (twk/two) header only / no header [null]\n"
	"  -O char  output type: b for TWO format, n for tab-delimited format (default: b)\n"
	"  -N       output in tab-delimited text format (see -O) [null]\n"
	"  -B       output in binary TWO/TWK format (see -O, default)[null]\n"
	"  -s       Hide all program messages [null]\n\n"

	// Twk parameters
	"Twk parameters\n"
	"  -G       drop genotypes in output [null]\n\n"

	// Two parameters
	"Two parameters\n"
	"  -r, --minR2  FLOAT   Pearson's R-squared minimum cut-off value (default: 0.1)\n"
	"  -R, --maxR2  FLOAT   Pearson's R-squared maximum cut-off value (default: 1.0)\n"
	"  -p, --minP   FLOAT   smallest P-value (default: 0)\n"
	"  -P, --maxP   FLOAT   largest P-value (default: 1)\n"
	"  -d, --minD   FLOAT   smallest D value (default: -1)\n"
	"  -D, --maxD   FLOAT   largest D value (default: 1)\n"
	"  -b, --minDP  FLOAT   smallest D' value (default: 0)\n"
	"  -B, --maxDP  FLOAT   largest D' value (default: 1)\n"
	"  -x, --minChi FLOAT   smallest Chi-squared CV (default: 0)\n"
	"  -X, --maxChi FLOAT   largest Chi-squared CV (default: inf)\n"
	"  -a, --minMHC FLOAT   minimum minor-haplotype count (default: 0)\n"
	"  -A, --maxMHC FLOAT   maximum minor-haplotype count (default: inf)\n"
	"  -m, --minMP  FLOAT   smallest model Chi-squared CV (default: 0)\n"
	"  -M, --maxMP  FLOAT   largest model Chi-squared CV (default: inf)\n"
	"  -f           INT     include FLAG value\n"
	"  -F           INT     exclude FLAG value\n"
	"  --min<cell>  FLOAT   smallest cell count (default: 0)\n"
	"  --max<cell>  FLOAT   largest cell count (default: inf)\n";
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
		{"minDP",	optional_argument, 0, 'd' },
		{"maxDP",	optional_argument, 0, 'D' },
		{"minChi",	optional_argument, 0, 'x' },
		{"maxChi",	optional_argument, 0, 'X' },
		{"minAlelles",	optional_argument, 0, 'a' },
		{"maxAlleles",	optional_argument, 0, 'A' },
		{"minMP",	optional_argument, 0, 'm' },
		{"maxMP",	optional_argument, 0, 'M' },
		{"flagInclude",	optional_argument, 0, 'f' },
		{"flagExclude",	optional_argument, 0, 'F' },
		{"headerOnly",	no_argument, 0, 'H' },
		{"noHeader",	no_argument, 0, 'h' },
		{"dropGenotypes",	optional_argument, 0, 'G' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	// Parameter defaults
	std::string input, output;
	Tomahawk::OutputFilter two_filter;
	bool outputHeader = true;
	int outputType = 1;
	bool dropGenotypes = false;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:o:r:R:p:P:d:D:x:X:a:A:m:M:f:F:HhGsBN", long_options, &long_index)) != -1){
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
			two_filter.minR2 = atof(optarg);
			if(two_filter.minR2 < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter r cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.minR2 > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter r has to be in range 0 < r < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'R':
			two_filter.maxR2 = atof(optarg);
			if(two_filter.maxR2 < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter R cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.maxR2 > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter R has to be in range 0 < R < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'd':
			two_filter.minDprime = atof(optarg);
			if(two_filter.minDprime < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter d cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.minDprime > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter d has to be in range 0 < d < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'D':
			two_filter.maxDprime = atof(optarg);
			if(two_filter.maxDprime < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter D cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.maxDprime > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter D has to be in range 0 < D < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'p':
			two_filter.minP = atof(optarg);
			if(two_filter.minP < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter p cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.minP > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter p has to be in range 0 < p < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'P':
			two_filter.maxP = atof(optarg);
			if(two_filter.maxP < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter P cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.maxP > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter P has to be in range 0 < P < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'x':
			two_filter.minChiSquared = atof(optarg);
			if(two_filter.minChiSquared < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter x cannot be negative" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'X':
			two_filter.maxChiSquared = atof(optarg);
			if(two_filter.maxChiSquared <= 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter X cannot be <= 0" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'm':
			two_filter.minPmodel = atof(optarg);
			if(two_filter.minPmodel < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter m cannot be negative" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'M':
			two_filter.maxPmodel = atof(optarg);
			if(two_filter.maxPmodel <= 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter M cannot be <= 0" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'a':
			two_filter.minMHF = atof(optarg);
			if(two_filter.minMHF < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter a cannot be negative" << std::endl;
				return(1);
			}
			break;

		case 'A':
			two_filter.maxMHF = atof(optarg);
			if(two_filter.maxMHF < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Parameter A cannot be negative" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'f':
			two_filter.filterValueInclude = atoi(optarg);
			two_filter.trigger();
			break;
		case 'F':
			two_filter.filterValueExclude = atoi(optarg);
			two_filter.trigger();
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

		//std::cerr << param << std::endl;

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
		reader.setWriteHeader(outputHeader);
		Tomahawk::OutputFilter& filter = reader.getFilter();
		filter = Tomahawk::OutputFilter(two_filter); // use copy ctor to transfer data

		if(!reader.setWriterType(outputType))
			return 1;

		if(!reader.Open(input))
			return 1;

		if(!reader.addRegions(filter_regions)){
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
