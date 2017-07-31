#include <getopt.h>

#include "utility.h"
#include "tomahawk/TomahawkReader.h"
#include "totempole/TotempoleReader.h"
#include "io/TomahawkOutput/TomahawkOutputReader.h"

void view_usage(void){
	programMessage();
	std::cout <<
	"Usage: " << Tomahawk::Constants::PROGRAM_NAME << " view [options] <in.twk>|<in.two>\n"
	"\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -o FILE  output file (- for stdout)\n\n"

	"Twk parameters\n"
	"  -G       (twk) drop genotypes in output [null]\n"
	"  -h/H     (twk/two) header only / no header [null]\n"
	"  -N       output in tab-delimited text format [null]\n"
	"  -B       output in binary TWO/TWK format (default)[null]\n"
	"  -t INT   number of CPU threads (default: maximum available)\n\n"

	// Two parameters
	"Two parameters\n"
	"  -p float number of parts to split problem into (default: 1)\n"
	"  -P float chosen part to compute (0 < -C < -c; default: 1)\n"
	"  -d float number of parts to split problem into (default: 1)\n"
	"  -D float chosen part to compute (0 < -C < -c; default: 1)\n"
	"  -a INT   minimum number of non-major genotypes in 2-by-2 matrix (default: 5)\n"
	"  -A INT   minimum number of non-major genotypes in 2-by-2 matrix (default: 5)\n"
	"  -f INT   minimum number of non-major genotypes in 2-by-2 matrix (default: 5)\n"
	"  -F INT   minimum number of non-major genotypes in 2-by-2 matrix (default: 5)\n"
	"  -r FLOAT Pearson's R-squared minimum cut-off value (default: 0.1)\n"
	"  -R FLOAT Pearson's R-squared maximum cut-off value (default: 1.0)\n"
	"  -d       Show real-time progress update in cerr [null]\n"
	"  -s       Hide all program messages [null]\n";
}

int view(int argc, char** argv){
	//argc -= 1; argv += 1;

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
	S32 minAlleles = 0, maxAlleles = std::numeric_limits<S32>::max();
	U16 flagInclude = 0, flagExclude = 0;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:o:P:p:a:A:R:r:f:F:d:D:w:W:hHGs", long_options, &long_index)) != -1){
		//std::cerr << c << ":" << (char)c << '\t' << long_index << std::endl;
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
	std::transform(end.begin(), end.end(), end.begin(), ::tolower);
	// std::cerr << end << std::endl;

	// Todo: action
	// Parse remainder parameters
	// Assume these parameters are contig or position values for filtering
	std::vector<std::string> filter_regions;
	for(U32 i = 2+hits; i < argc; ++i){
		std::string param(&argv[i][0]);
		//std::cerr << i << '\t' << param << std::endl;

		if(!Tomahawk::Helpers::parsePositionalStringTWO(param)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Illegal parse of input string: " << param << std::endl;
			return 1;
		}

		filter_regions.push_back(param);
	}

	if(end == Tomahawk::Constants::OUTPUT_SUFFIX){
		Tomahawk::TomahawkReader tomahawk;
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
		if(!filter.setFilterRsquared(minR2, maxR2)) return false;
		filter.setFilterInclude(flagInclude);
		filter.setFilterExclude(flagExclude);
		//if(!filter.setFilterMGF(minAlleles)) return false;
		if(!filter.setFilterP(minP, maxP)) return false;
		if(!filter.setFilterDprime(minDprime, maxDprime)) return false;

		if(!reader.Open(input)){
			return false;
		}

		reader.AddRegions(filter_regions);

		if(!reader.view(input)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Failed to read!" << std::endl;
			return 1;
		}
	} else {
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Unrecognized input file format: " << input << std::endl;
		return 1;
	}

	return 0;
}
