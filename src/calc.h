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
#include <regex>

#include "support/MagicConstants.h"
#include "tomahawk/TomahawkCalc.h"

void calc_usage(void){
	programMessage();
	std::cerr <<
	"About:  Calculate linkage disequilibrium\n"
	"        Force phased -p or unphased -u for faster calculations if\n"
	"        the all variant sites are guaranteed to have the given phasing.\n"
	"Usage:  " << tomahawk::constants::PROGRAM_NAME << " calc [options] -i <in.twk> -o <output.two>\n\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -o FILE  output file (required)\n"
	"  -t INT   number of CPU threads (default: maximum available)\n"
	"  -c INT   number of parts to split problem into (must be in c!2 + c unless -w is triggered)\n"
	"  -C INT   chosen part to compute (0 < -C < -c)\n"
	"  -w INT   sliding window width in bases (approximative)\n"
	"  -p       force computations to use phased math [null]\n"
	"  -u       force computations to use unphased math [null]\n"
	"  -f       use fast-mode: output will have correlations only (no matrices or tests) [null]\n"
	"  -S INT   number of individuals to sample in fast-mode (default: 1000)\n"
	"  -a INT   minimum number of non-major genotypes in 2-by-2 matrix (default: 1)\n"
	"  -P FLOAT Fisher's exact test / Chi-squared cutoff P-value (default: 1)\n"
	"  -r FLOAT Pearson's R-squared minimum cut-off value (default: 0.1)\n"
	"  -R FLOAT Pearson's R-squared maximum cut-off value (default: 1.0)\n"
	"  -d       Show real-time progress update in cerr [null]\n"
	"  -s       Hide all program messages [null]\n";
}

int calc(int argc, char** argv){
	//argc -= 2; argv += 2;

	int c;

	if(argc < 3){
		calc_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",             required_argument, 0, 'i' },
		{"threads",           optional_argument, 0, 't' },
		{"output",            required_argument, 0, 'o' },
		{"parts",             optional_argument, 0, 'c' },
		{"partStart",         optional_argument, 0, 'C' },
		{"minP",              optional_argument, 0, 'P' },
		{"phased",            no_argument,       0, 'p' },
		{"unphased",          no_argument,       0, 'u' },
		{"fast-mode",         no_argument,       0, 'f' },
		{"samples",           optional_argument, 0, 'S' },
		{"minR2",             optional_argument, 0, 'r' },
		{"maxR2",             optional_argument, 0, 'R' },
		{"minMHF",            optional_argument, 0, 'a' },
		{"maxMHF",            optional_argument, 0, 'A' },
		{"detailedProgress",  no_argument,       0, 'd' },
		{"silent",            no_argument,       0, 's' },
		// Not implemented
		{"windowBases",       optional_argument, 0, 'w' },
		{"windowSites",       optional_argument, 0, 'W' },
		{"longHelp",          optional_argument, 0, '?' },
		{0,0,0,0}
	};

	tomahawk::TomahawkCalc tomahawk;
	tomahawk::TomahawkCalcParameters& parameters = tomahawk.getParameters();
	std::string input;
	std::string output;

	double windowBases = -1, windowPosition = -1; // not implemented

	while ((c = getopt_long(argc, argv, "i:o:t:puP:a:A:r:R:w:W:S:sdc:C:f?", long_options, &option_index)) != -1){
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
		case 't':
			parameters.n_threads = atoi(optarg);
			if(parameters.n_threads <= 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a non-positive number of worker threads" << std::endl;
				return(1);
			}
			break;
		case 'c':
			parameters.n_chunks = atoi(optarg);
			if(parameters.n_chunks <= 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a negative or zero amount of partitions" << std::endl;
				return(1);
			}
			break;
		case 'C':
			parameters.chunk_selected = atoi(optarg);
			--parameters.chunk_selected;
			if(parameters.chunk_selected < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a non-positive start partition" << std::endl;
				return(1);
			}
			break;
		case 'r':
			parameters.R2_min = atof(optarg);
			if(parameters.R2_min < 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a negative minimum R-squared value" << std::endl;
				return(1);
			} else if(parameters.R2_min > 1){
				std::cerr << tomahawk::helpers::timestamp("ERROR")<< "Cannot have minimum R-squared value > 1" << std::endl;
				return(1);
			}
			break;
		case 'R':
		parameters.R2_max = atof(optarg);
		if(parameters.R2_max < 0){
			std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a negative maximum R-squared value" << std::endl;
		return(1);
		} else if(parameters.R2_max > 1){
			std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a maximum R-squared value > 1" << std::endl;
		return(1);
		}
		break;

		case 'p':
		  parameters.force = tomahawk::TomahawkCalcParameters::force_method::phasedFunction;
		  break;

		case 'u':
		  parameters.force = tomahawk::TomahawkCalcParameters::force_method::unphasedFunction;
		  break;

		case 'f':
		  parameters.fast_mode = true;
		  break;

		case 'P':
		  parameters.P_threshold = atof(optarg);
		  if(parameters.P_threshold < 0){
			  std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a negative cutoff P-value" << std::endl;
			return(1);
		  } else if(parameters.P_threshold > 1){
			  std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a cutoff P-value > 1" << std::endl;
			return(1);
		  }
		  break;
		case 'a':
		parameters.minimum_sum_alternative_haplotype_count = atoi(optarg);
		if(parameters.minimum_sum_alternative_haplotype_count < 0){
			std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have negative minimum allele count" << std::endl;
			return(1);
		}
		break;

		case 'A':
		parameters.maximum_sum_alternative_haplotype_count = atoi(optarg);
		if(parameters.maximum_sum_alternative_haplotype_count < 0){
			std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have negative maximum allele count" << std::endl;
			return(1);
		}
		break;

		case 'w':
		  if(std::regex_match(std::string(optarg), std::regex("^(([0-9]+)|([0-9]+[eE]{1}[0-9]+))$")) == false){
			  std::cerr << "not an integer" << std::endl;
			  return(1);
		  }
		  if(std::regex_match(std::string(optarg), std::regex("^[0-9]+$"))){
			windowBases = atoi(optarg);
		  } else if(std::regex_match(std::string(optarg), std::regex("^[0-9]+[eE]{1}[0-9]+$"))){
			  windowBases = atof(optarg);
		  }

		  if(windowBases <= 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a non-positive window size" << std::endl;
				return(1);
			}
		parameters.window_mode = true;
		parameters.n_window_bases = windowBases;
		break;

		case 'W':
			windowPosition = atoi(optarg);
			if(windowPosition <= 0){
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Cannot have a non-positive window size" << std::endl;
				return(1);
			}

		break;

		case 's':
			SILENT = 1;
			parameters.detailed_progress = false;
			break;

		case 'd':
		  SILENT = 0;
		  parameters.detailed_progress = true;
		  break;

		default:
		  std::cerr << tomahawk::helpers::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
		  return(1);
		}
	}

	if(input.length() == 0){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(output.length() == 0){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "No output value specified..." << std::endl;
		return(1);
	}

	// Print messages
	if(!SILENT){
		programMessage();
		std::cerr << tomahawk::helpers::timestamp("LOG") << "Calling calc..." << std::endl;
	}


	// Parse Tomahawk
	if(!tomahawk.Open(input, output)){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "Failed build!" << std::endl;
		return 1;
	}

	return(tomahawk.Calculate());
}
