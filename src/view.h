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

#include "tomahawk/tomahawk_reader.h"
#include "tomahawk/output_filter.h"
#include "tomahawk/tomahawk_output_reader.h"
#include "utility.h"

void view_usage(void){
	programMessage();
	std::cerr <<
	"About:  Convert binary TWK->VCF or TWO->LD, subset and slice TWK/TWO data\n\n"
	"Usage:  " << tomahawk::constants::PROGRAM_NAME << " view [options] -i <in.two>\n\n"
	"Options:\n"
	"  -i FILE   input Tomahawk (required)\n"
	"  -h/H      (twk/two) header only / no header\n"
	"  -N        output in tab-delimited text format (see -O)\n"
	"  -B        output in binary TWO/TWK format (see -O, default)\n"
	"  -I STRING filter interval <contig>:pos-pos (TWK/TWO) or linked interval <contig>:pos-pos,<contig>:pos-pos (TWO only)\n"
	//"  -J        output JSON object\n"
	"  -s        Hide all program messages\n\n"

	// Twk parameters
	"TWK parameters\n"
	"  -G       drop genotypes in output\n\n"

	// Two parameters
	"TWO parameters\n"
	"  -o FILE   output file (- for stdout; default: -)\n"
	"  -O char   output type: b for TWO format, u for tab-delimited LD format\n"
	"  -r, --minR2  FLOAT   Pearson's R-squared minimum cut-off value\n"
	"  -R, --maxR2  FLOAT   Pearson's R-squared maximum cut-off value\n"
	"  -z, --minR   FLOAT   Pearson's R minimum cut-off value\n"
	"  -Z, --maxR   FLOAT   Pearson's R maximum cut-off value\n"
	"  -p, --minP   FLOAT   smallest P-value (default: 0)\n"
	"  -P, --maxP   FLOAT   largest P-value (default: 1)\n"
	"  -d, --minD   FLOAT   smallest D value (default: -1)\n"
	"  -D, --maxD   FLOAT   largest D value (default: 1)\n"
	"  -b, --minDP  FLOAT   smallest D' value (default: 0)\n"
	"  -B, --maxDP  FLOAT   largest D' value (default: 1)\n"
	"  -x, --minChi FLOAT   smallest Chi-squared CV of contingency table (default: 0)\n"
	"  -X, --maxChi FLOAT   largest Chi-squared CV of contingency table (default: inf)\n"
	"  -a, --minMHC FLOAT   minimum minor-haplotype count (default: 0)\n"
	"  -A, --maxMHC FLOAT   maximum minor-haplotype count (default: inf)\n"
	"  -m, --minMCV FLOAT   smallest Chi-squared CV of unphased model (default: 0)\n"
	"  -M, --maxMCV FLOAT   largest Chi-squared CV of unphased model (default: inf)\n"
	"  -1, --maxP1  FLOAT   maximum REF_REF count\n"
	"  -2, --maxP2  FLOAT   maximum REF_ALT count\n"
	"  -3, --maxQ1  FLOAT   maximum ALT_REF count\n"
	"  -4, --maxQ2  FLOAT   maximum ALT_ALT count\n"
	"  -f           INT     include FLAG value\n"
	"  -F           INT     exclude FLAG value\n"
	"  -u                   output only the upper triangular values\n"
	"  -l                   output only the lower triangular values\n";
}

int view(int argc, char** argv){
	if(argc < 3){
		view_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{"output",      optional_argument, 0, 'o' },
		{"output-type", optional_argument, 0, 'O' },
		{"minP",        optional_argument, 0, 'p' },
		{"maxP",        optional_argument, 0, 'P' },
		{"minR",        optional_argument, 0, 'z' },
		{"maxR",        optional_argument, 0, 'Z' },
		{"minR2",       optional_argument, 0, 'r' },
		{"maxR2",       optional_argument, 0, 'R' },
		{"minDP",       optional_argument, 0, 'b' },
		{"maxDP",       optional_argument, 0, 'B' },
		{"minD",        optional_argument, 0, 'd' },
		{"maxD",        optional_argument, 0, 'D' },
		{"minChi",      optional_argument, 0, 'x' },
		{"maxChi",      optional_argument, 0, 'X' },
		{"minAlelles",  optional_argument, 0, 'a' },
		{"maxAlleles",  optional_argument, 0, 'A' },
		{"minMCV",      optional_argument, 0, 'm' },
		{"maxMCV",      optional_argument, 0, 'M' },
		{"JSON",        optional_argument, 0, 'J' },
		{"flagInclude", optional_argument, 0, 'f' },
		{"flagExclude", optional_argument, 0, 'F' },
		{"upperTriangular", no_argument, 0, 'u' },
		{"lowerTriangular", no_argument, 0, 'l' },
		{"headerOnly",  no_argument, 0, 'H' },
		{"noHeader",    no_argument, 0, 'h' },
		{"dropGenotypes",optional_argument, 0, 'G' },
		{"interval",    optional_argument, 0, 'I' },
		{"silent",      no_argument, 0,  's' },
		{"maxP1",  optional_argument, 0, '1' },
		{"maxP2",  optional_argument, 0, '2' },
		{"maxQ1",  optional_argument, 0, '3' },
		{"maxQ2",  optional_argument, 0, '4' },

		{0,0,0,0}
	};

	// Parameter defaults
	std::string input, output;
	tomahawk::OutputFilter two_filter;
	tomahawk::TomahawkOutputReaderParameters two_parameters;
	bool outputHeader     = true;
	bool outputHeaderOnly = false;
	bool dropGenotypes    = false;
	std::vector<std::string> filter_regions;
	//bool output_JSON = false;
	std::string temp;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:o:O:r:R:p:P:d:D:x:X:a:A:m:M:f:F:I:HhGsb:B:Nulz:Z:1:2:3:4:", long_options, &long_index)) != -1){
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
			two_parameters.input_file = input;
			break;
		case 'o':
			output = std::string(optarg);
			two_parameters.output_file = output;
			break;
		case 'O':
			temp = std::string(optarg);
			if(temp.size() == 1 && temp[0] == 'b')      two_parameters.output_type = tomahawk::TWK_OUTPUT_TWO;
			else if(temp.size() == 1 && temp[0] == 'u') two_parameters.output_type = tomahawk::TWK_OUTPUT_LD;
			else {
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Unknown value " << temp << " for parameter -O..." << std::endl;
				return(1);
			}
			break;
		case 'I':
			filter_regions.push_back(std::string(optarg));
			break;
		//case 'J':
		//	output_JSON = true;
		//	--hits;
		//	break;
		case 'u':
			two_filter.setFilterUpperTriangular(true);
			two_filter.trigger();
			--hits;
			break;
		case 'l':
			two_filter.setFilterLowerTriangular(true);
			two_filter.trigger();
			--hits;
			break;
		case 'r':
			two_filter.minR2 = atof(optarg);
			if(two_filter.minR2 < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter r cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.minR2 > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter r has to be in range 0 < r < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'R':
			two_filter.maxR2 = atof(optarg);
			if(two_filter.maxR2 < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter R cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.maxR2 > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter R has to be in range 0 < R < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'z':
			two_filter.minR = atof(optarg);
			if(two_filter.minR < -1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter z cannot be < -1" << std::endl;
				return(1);
			}
			if(two_filter.minR > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter z has to be in range -1 < r < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'Z':
			two_filter.maxR = atof(optarg);
			if(two_filter.maxR < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter Z cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.maxR > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter Z has to be in range 0 < R < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'd':
			two_filter.minD = atof(optarg);
			if(two_filter.minD < -1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter d cannot be < -1" << std::endl;
				return(1);
			}
			if(two_filter.minD > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter d has to be in range -1 < d < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'D':
			two_filter.maxD = atof(optarg);
			if(two_filter.maxD < -1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter D cannot be < -1" << std::endl;
				return(1);
			}
			if(two_filter.maxD > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter D has to be in range -1 < D < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'p':
			two_filter.minP = atof(optarg);
			if(two_filter.minP < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter p cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.minP > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter p has to be in range 0 < p < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'P':
			two_filter.maxP = atof(optarg);
			if(two_filter.maxP < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter P cannot be negative" << std::endl;
				return(1);
			}
			if(two_filter.maxP > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter P has to be in range 0 < P < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'x':
			two_filter.minChiSquaredTable = atof(optarg);
			if(two_filter.minChiSquaredTable < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter x cannot be negative" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'X':
			two_filter.maxChiSquaredTable = atof(optarg);
			if(two_filter.maxChiSquaredTable <= 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter X cannot be <= 0" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'm':
			two_filter.minChiSquaredModel = atof(optarg);
			if(two_filter.minChiSquaredModel < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter m cannot be negative" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'M':
			two_filter.maxChiSquaredModel = atof(optarg);
			if(two_filter.maxChiSquaredModel <= 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter M cannot be <= 0" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'a':
			two_filter.minMHF = atof(optarg);
			if(two_filter.minMHF < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter a cannot be negative" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'A':
			two_filter.maxMHF = atof(optarg);
			if(two_filter.maxMHF < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter A cannot be negative" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'f':
			two_filter.FLAGInclude = atoi(optarg);
			two_filter.trigger();
			break;
		case 'F':
			two_filter.FLAGExclude = atoi(optarg);
			two_filter.trigger();
			break;
		case 's':
			SILENT = 1;
			--hits;
			break;
		case 'h':
			outputHeaderOnly = true;
			--hits;
			break;
		case 'H':
			outputHeader = false;
			--hits;
			break;

		case 'b':
			two_filter.minDprime = atof(optarg);
			if(two_filter.minDprime < -1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter b cannot be < -1" << std::endl;
				return(1);
			}
			if(two_filter.minDprime > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter b has to be in range -1 < d < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		case 'B':
			two_filter.maxDprime = atof(optarg);
			if(two_filter.maxDprime < -1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter B cannot be < -1" << std::endl;
				return(1);
			}
			if(two_filter.maxDprime > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter B has to be in range -1 < D < 1" << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case 'G':
			dropGenotypes = true;
			break;


		case '1':
			two_filter.maxP1 = atof(optarg);
			if(two_filter.maxP1 < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter 1 has to be positive..." << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case '2':
			two_filter.maxP2 = atof(optarg);
			if(two_filter.maxP2 < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter 2 has to be positive..." << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case '3':
			two_filter.maxQ1 = atof(optarg);
			if(two_filter.maxQ1 < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter 3 has to be positive..." << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;

		case '4':
			two_filter.maxQ2 = atof(optarg);
			if(two_filter.maxQ2 < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Parameter 4 has to be positive..." << std::endl;
				return(1);
			}
			two_filter.trigger();
			break;
		}
	}

	if(input.length() == 0){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "No input file specified..." << std::endl;
		return(1);
	}

	if(!SILENT){
		programMessage();
		std::cerr << tomahawk::helpers::timestamp("LOG") << "Calling view..." << std::endl;
	}

	std::vector<std::string> inputFile_parts = tomahawk::helpers::split(input, '.');
	std::string& end = inputFile_parts[inputFile_parts.size() - 1];
	std::transform(end.begin(), end.end(), end.begin(), ::tolower); // transform chars to lower case

	if(end == tomahawk::constants::OUTPUT_SUFFIX){
		tomahawk::TomahawkReader tomahawk;
		tomahawk.setDropGenotypes(dropGenotypes);
		tomahawk.setShowHeader(outputHeader);
		if(!tomahawk.open(input)){
			std::cerr << tomahawk::helpers::timestamp("ERROR") << "Failed to open!" << std::endl;
			return 1;
		}

		if(outputHeaderOnly){
			tomahawk.printHeader(std::cout);
			return(0);
		}

		if(!tomahawk.addRegions(filter_regions)){
			std::cerr << tomahawk::helpers::timestamp("ERROR") << "Failed to add region!" << std::endl;
			return 1;
		}

		//tomahawk.summaryIndividuals();
		tomahawk.outputBlocks();

	} else if(end == tomahawk::constants::OUTPUT_LD_SUFFIX){
		tomahawk::TomahawkOutputReader reader;
		reader.parameters_ = two_parameters;
		tomahawk::OutputFilter& filter = reader.getFilter();
		filter = tomahawk::OutputFilter(two_filter); // use copy ctor to transfer data

		if(!reader.open(input))
			return 1;

		if(outputHeaderOnly){
			reader.printHeader(std::cout);
			return(0);
		}

		if(!reader.addRegions(filter_regions)){
			std::cerr << tomahawk::helpers::timestamp("ERROR") << "Failed to add region!" << std::endl;
			return 1;
		}

		if(!reader.view())
			return 1;

	} else {
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "Unrecognised input file format: " << input << std::endl;
		return 1;
	}

	return 0;
}
