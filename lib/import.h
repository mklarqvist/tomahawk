#ifndef IMPORT_H_
#define IMPORT_H_

#include "getopt.h"
#include "utility.h"
#include "tomahawk.h"
#include "importer.h"

void import_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Convert BCF->TWK/; subset and slice TWK/TWO data\n"
	"        Only biallelic diploid genotypes from SNVs will be retained\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " import [options] -i <in.bcf> -o <output.twk>\n\n"
	"Options:\n"
	"  -i FILE  input BCF file (required)\n"
	"  -o FILE  output file prefix (required)\n"
	"  -f       Flip reference and alternative alleles when major is the alternative allele\n"
	"  -r       Do NOT filter out variant sites that are univariate for REF or ALT\n"
	"  -n FLOAT Missingness fraction cutoff (default: 0.95)\n"
	"  -b INT   Block size (default: 500)\n"
	"  -L INT   Compression level in range 1-20 (default: 1)\n"
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
		{"input",       required_argument, 0,  'i' },
		{"output",      optional_argument, 0,  'o' },
		//{"extend",      optional_argument, 0,  'e' },
		{"filter-univariate", optional_argument, 0,  'r' },
		{"flip",        optional_argument, 0,  'f' },
		{"missingness", optional_argument, 0,  'n' },
		{"compression-level", optional_argument, 0,  'L' },
		{"block-size", optional_argument, 0,  'b' },
		{"hwe", optional_argument, 0,  'H' },
		{0,0,0,0}
	};
	tomahawk::twk_vimport_settings settings;

	while ((c = getopt_long(argc, argv, "i:o:rfn:b:L:H:?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			settings.input = std::string(optarg);
			break;
		case 'o':
			settings.output = std::string(optarg);
			break;
		//case 'e':
		//	extension_mode = true;
		//	extend = std::string(optarg);
		//	break;
		case 'n':
			settings.threshold_miss = atof(optarg);
			if(settings.threshold_miss < 0){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot set missingness filter to < 0..." << std::endl;
				return(1);
			}

			if(settings.threshold_miss > 1){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot set missingness filter to > 1..." << std::endl;
				return(1);
			}

			break;
		case 'H':
			settings.hwe = atof(optarg);
			if(settings.hwe < 0){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot set Hardy-Weinberg filter to < 0..." << std::endl;
				return(1);
			}

			if(settings.hwe > 1){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot set Hardy-Weinberg filter to > 1..." << std::endl;
				return(1);
			}

			break;
		case 'r':
			settings.remove_univariate = false;
			break;
		case 'f':
			settings.flip_major_minor = false;
			break;
		case 'b':
			settings.block_size = atoi(optarg);
			break;
		case 'L':
			settings.c_level = atoi(optarg);
			break;

		default:
			std::cerr << tomahawk::utility::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if(settings.input.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(settings.output.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No output value specified..." << std::endl;
		return(1);
	}

	// Print messages
	tomahawk::ProgramMessage();
	std::cerr << tomahawk::utility::timestamp("LOG") << "Calling import..." << std::endl;

	tomahawk::twk_variant_importer importer;
	if(importer.Import(settings) == false){
		std::cerr << "failed import" << std::endl;
		return 1;
	}

	return 0;
}

#endif /* IMPORT_H_ */
