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
#include <regex>

#include "ld.h"

void calc_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Calculate linkage disequilibrium\n"
	"        Force phased -p or unphased -u for faster calculations if\n"
	"        all variant sites are guaranteed to have the given phasing.\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " calc [options] -i <in.twk> -o <output.two>\n\n"
	"Options:\n"
	"  -i FILE   input Tomahawk (required)\n"
	"  -o FILE   output file or file prefix (required)\n"
	"  -t INT    number of CPU threads (default: maximum available)\n"
	"  -c INT    number of subproblems to split compute into (must be in (c!2 + c))\n"
	"  -C INT    chosen part to compute (0 < -C < -c)\n"
	"  -m        run in low-memory mode: this is considerably slower but use no more memory than\n"
	"               block1*variants + block2*variants\n"
	"  -M        use phased bitmaps in low-memory mode. Automatically triggers -m and -p.\n"
	"  -b        number of records in a block. Has an effect on memory usage only when -m is set.\n"
	"  -w INT    sliding window width in bases\n"
	"  -I STRING filter interval <contig>:pos-pos (see manual)\n"
	"  -p        force computations to use phased math\n"
	"  -u        force computations to use unphased math\n"
	//"  -S INT    trigger sampling mode: number of individuals to sample when allele counts are large.\n"
	"  -P FLOAT  Fisher's exact test / Chi-squared cutoff P-value (default: 1)\n"
	"  -r FLOAT  Pearson's R-squared minimum cut-off value (default: 0.1)\n"
	//"  -R FLOAT  Pearson's R-squared maximum cut-off value (default: 1.0)\n"
	"  -k INT    compression level to use (default: 1, max = 22).\n" << std::endl;
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
		{"interval",          optional_argument, 0, 'I' },
		{"parts",             optional_argument, 0, 'c' },
		{"partStart",         optional_argument, 0, 'C' },
		{"low-memory",        optional_argument, 0, 'm' },
		{"block-size",        optional_argument, 0, 'b' },
		{"bitmaps",           optional_argument, 0, 'M' },
		{"compression-level", optional_argument, 0, 'k' },

		{"cross-chr-only",    no_argument, 0, 'X' },
		{"no-cross-chr",      no_argument, 0, 'x' },

		{"minP",              optional_argument, 0, 'P' },
		{"force-phased",      no_argument,       0, 'p' },
		{"force-unphased",    no_argument,       0, 'u' },
		{"samples",           optional_argument, 0, 'S' },
		{"minR2",             optional_argument, 0, 'r' },
		//{"maxR2",             optional_argument, 0, 'R' },

		{"detailedProgress",  no_argument,       0, 'd' },
		{"silent",            no_argument,       0, 's' },
		// Not implemented
		{"windowBases",       optional_argument, 0, 'w' },
		{0,0,0,0}
	};

	tomahawk::twk_ld_settings settings;
	//std::vector<std::string> filter_regions;

	while ((c = getopt_long(argc, argv, "i:o:t:puP:a:A:r:w:S:I:sdc:C:mMb:xXk:?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			settings.in = std::string(optarg);
			break;
		case 'o':
			settings.out = std::string(optarg);
			break;
		case 'I':
			settings.ival_strings.push_back(std::string(optarg));
			break;
		case 'm':
			settings.low_memory = true;
			break;
		case 'p':
			settings.force_phased = true;
			settings.forced_unphased = false;
			break;
		case 'u':
			settings.forced_unphased = true;
			settings.force_phased = false;
			break;
		case 'M':
			settings.force_phased = true;
			settings.low_memory = true;
			settings.bitmaps = true;
			break;
		case 't':
			settings.n_threads = atoi(optarg);
			if(settings.n_threads <= 0){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a non-positive number of worker threads" << std::endl;
				return(1);
			}
			break;

		case 'b':
			settings.bl_size = atoi(optarg);
			if(settings.bl_size <= 0){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a non-positive number of entries in a block!" << std::endl;
				return(1);
			}
			break;
		case 'c':
			settings.n_chunks = atoi(optarg);
			if(settings.n_chunks <= 0){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a negative or zero amount of partitions" << std::endl;
				return(1);
			}
			break;
		case 'C':
			settings.c_chunk = atoi(optarg);
			--settings.c_chunk;
			if(settings.c_chunk < 0){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a non-positive start partition" << std::endl;
				return(1);
			}
			break;
		case 'r':
			settings.minR2 = atof(optarg);
			if(settings.minR2 < 0){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a negative minimum R-squared value" << std::endl;
				return(1);
			} else if(settings.minR2 > 1){
				std::cerr << tomahawk::utility::timestamp("ERROR")<< "Cannot have minimum R-squared value > 1" << std::endl;
				return(1);
			}
			break;
		/*
		case 'R':
			settings.maxR2 = atof(optarg);
			if(settings.maxR2 < 0){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a negative maximum R-squared value" << std::endl;
			return(1);
			} else if(settings.maxR2 > 1){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a maximum R-squared value > 1" << std::endl;
			return(1);
			}
			break;
		*/
		case 'P':
		  settings.minP = atof(optarg);
		  if(settings.minP < 0){
			  std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a negative cutoff P-value" << std::endl;
			return(1);
		  } else if(settings.minP > 1){
			  std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a cutoff P-value > 1" << std::endl;
			return(1);
		  }
		  break;

		case 'w':
			settings.window = true;
		  if(std::regex_match(std::string(optarg), std::regex("^(([0-9]+)|([0-9]+[eE]{1}[0-9]+))$")) == false){
			  std::cerr << "not an integer" << std::endl;
			  return(1);
		  }
		  if(std::regex_match(std::string(optarg), std::regex("^[0-9]+$"))){
			settings.l_window = atoi(optarg);
		  } else if(std::regex_match(std::string(optarg), std::regex("^[0-9]+[eE]{1}[0-9]+$"))){
			  settings.l_window = atof(optarg);
		  }

		  if( settings.l_window <= 0){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a non-positive window size" << std::endl;
				return(1);
			}

		  break;

		case 'k':
		settings.c_level = std::atoi(optarg);
		break;


		default:
		  std::cerr << tomahawk::utility::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
		  return(1);
		}
	}

	if(settings.in.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(settings.out.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No output value specified..." << std::endl;
		return(1);
	}

	// Print messages
	tomahawk::ProgramMessage();
	std::cerr << tomahawk::utility::timestamp("LOG") << "Calling calc..." << std::endl;


	tomahawk::twk_ld ld;
	if(ld.Compute(settings) == false) return 1;
	return 0;
}
