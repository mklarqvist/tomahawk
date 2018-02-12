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
#include "tomahawk/TomahawkImporter.h"
#include "utility.h"

void import_usage(void){
	programMessage();
	std::cerr <<
	"About:  Convert VCF/BCF->TWK/; subset and slice TWK/TWO data\n"
	"        Only biallelic diploid genotypes from SNVs will be retained\n"
	"Usage:  " << Tomahawk::Constants::PROGRAM_NAME << " import [options] -i <in.vcf>/<in.bcf> -o <output.twk>\n\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -o FILE  output file prefix (required)\n"
	"  -h FLOAT Hardy-Weinberg P-value cutoff (default: 0)\n"
	"  -m FLOAT Minor-allele frequency (MAF) cutoff (default: 0)\n"
	"  -n FLOAT Missingness percentage cutoff (default: 0.2)\n"
	"  -s       Hide all program messages [null]\n";
}

int import(int argc, char** argv){
	int c;
	if(argc < 3){
		import_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",		required_argument, 0,  'i' },
		{"output",		optional_argument, 0,  'o' },
		{"extend",		optional_argument, 0,  'e' },
		{"hwep",		optional_argument, 0,  'h' },
		{"missingness",		optional_argument, 0,  'n' },
		{"maf",		optional_argument, 0,  'm' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	std::string input;
	std::string output;
	std::string extend;
	bool extension_mode = false;
	SILENT = 0;
	double hwe_p = 0;
	double maf = 0;
	double missingness = 0.2;

	while ((c = getopt_long(argc, argv, "i:o:e:h:m:n:s?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			input = std::string(optarg);
			break;
		case 'o':
			output = std::string(optarg);
			break;
		case 'e':
			extension_mode = true;
			extend = std::string(optarg);
			break;
		case 'h':
			hwe_p = atof(optarg);
			if(hwe_p < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Cannot set Hardy-Weinberg filter to < 0..." << std::endl;
				return(1);
			}

			if(hwe_p > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Cannot set Hardy-Weinberg filter to > 1..." << std::endl;
				return(1);
			}

			break;
		case 'n':
			missingness = atof(optarg);
			if(missingness < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Cannot set missingness filter to < 0..." << std::endl;
				return(1);
			}

			if(missingness > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Cannot set missingness filter to > 1..." << std::endl;
				return(1);
			}

			break;
		case 'm':
			maf = atof(optarg);
			if(maf < 0){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Cannot set MAF filter to < 0..." << std::endl;
				return(1);
			}

			if(maf > 1){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Cannot set MAF filter to > 1..." << std::endl;
				return(1);
			}

			break;
		case 's':
			SILENT = 1;
			break;

		default:
			std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if(input.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(!extension_mode && output.length() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No output value specified..." << std::endl;
		return(1);
	}

	if(extension_mode && extend.size() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "No file to extend provided..." << std::endl;
		return(1);
	}

	// Print messages
	if(!SILENT){
		programMessage();
		std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling import..." << std::endl;
	}


	Tomahawk::TomahawkImporter importer(input, output);
	importer.getFilters().HWE_P = hwe_p;
	importer.getFilters().MAF = maf;
	importer.getFilters().missingness = missingness;

	return true;

	if(!extension_mode){
		if(!importer.Build())
			return 1;

	} else {
		if(!importer.Extend(extend))
			return 1;

	}

	return 0;
}
