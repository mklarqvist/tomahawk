#include <getopt.h>

#include "utility.h"
#include "tomahawk/MagicConstants.h"
#include "tomahawk/TotempoleReader.h"
#include "tomahawk/TomahawkReader.h"

#define CALC_DEFAULT_MINR2 0.1
#define CALC_DEFAULT_MAXR2 1
#define CALC_DEFAULT_MINP  1e-4
#define CALC_DEFAULT_MINALLELES 5
#define CALC_DEFAULT_MAXALLELES std::numeric_limits<U64>::max()

void calc_usage(void){
	programMessage();
	std::cout <<
	"Usage: " << Tomahawk::Constants::PROGRAM_NAME << " calc [options] <in.twk>\n"
	"\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -o FILE  output file (- for stdout)\n"
	"  -N       output in tab-delimited text format [null]\n"
	"  -B       output in binary TWO format (default)[null]\n"
	"  -t INT   number of CPU threads (default: maximum available)\n"
	"  -c INT   number of parts to split problem into (default: 1)\n"
	"  -C INT   chosen part to compute (0 < -C < -c; default: 1)\n"
	"  -p       force computations to use phased math [null]\n"
	"  -u       force computations to use unphased math [null]\n"
	"  -a INT   minimum number of non-major genotypes in 2-by-2 matrix (default: 5)\n"
	"  -P FLOAT Fisher's exact test / Chi-squared cutoff P-value (default: 1e-4)\n"
	"  -r FLOAT Pearson's R-squared minimum cut-off value (default: 0.1)\n"
	"  -R FLOAT Pearson's R-squared maximum cut-off value (default: 1.0)\n"
	"  -S       Hide all program messages [null]\n";
}

int calc(int argc, char** argv){
	typedef Tomahawk::IO::TomahawkCalculationWriterInterace::compression compression_type;

	//argc -= 2; argv += 2;

	int c;

	if(argc < 3){
		calc_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",		required_argument, 0,  'i' },
		{"threads",		optional_argument,       0,  't' },
		{"output",		required_argument, 0,  'o' },
		{"parts",		optional_argument, 0,  'c' },
		{"partStart",	optional_argument, 0,  'C' },
		{"natural",		no_argument, 0,  'N' },
		{"binary",		no_argument, 0,  'B' },
		{"minP",		optional_argument, 0,  'P' },
		{"phased",		no_argument, 0,  'p' },
		{"unphased",	no_argument, 0,  'u' },
		{"minR2",		optional_argument, 0,  'r' },
		{"maxR2",		optional_argument, 0,  'R' },
		{"minalelles",	optional_argument, 0,  'a' },
		{"maxalleles",	optional_argument, 0,  'A' },
		{"silent",		no_argument, 0,  's' },
		// Not implemented
		{"windowBases",	optional_argument, 0,  'w' },
		{"windowPosition",optional_argument, 0,  'W' },
		{"longHelp",	optional_argument, 0, '?' },
		{0,0,0,0}
	};

	std::string input;
	std::string output = "-";
	compression_type type = compression_type::binary;
	S32 threads = std::thread::hardware_concurrency();
	S32 parts = 1;
	S32 startPart = 0;
	float minR2 = CALC_DEFAULT_MINR2;
	float maxR2 = CALC_DEFAULT_MAXR2;
	double minP = CALC_DEFAULT_MINP;
	int64_t minAlleles = CALC_DEFAULT_MINALLELES;
	int64_t maxAlleles = CALC_DEFAULT_MAXALLELES;
	S32 windowBases = -1, windowPosition = -1; // not implemented
	bool silent = false;
	bool phased = false;
	bool forceFunction = false;

	while ((c = getopt_long(argc, argv, "i:o:t:puP:a:A:r:R:w:W:sNBc:C:?", long_options, &option_index)) != -1){
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
			threads = atoi(optarg);
			if(threads <= 0){
				std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a non-positive number of worker threads" << std::endl;
				return(1);
			}
			break;
		case 'N':
			type = compression_type::natural;
			break;
		case 'B':
			type = compression_type::binary;
			break;
		case 'c':
			parts = atoi(optarg);
			if(parts <= 0){
				std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a negative or zero amount of partitions" << std::endl;
				return(1);
			}
			break;
		case 'C':
			startPart = atoi(optarg);
			--startPart;
			if(startPart < 0){
				std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a non-positive start partition" << std::endl;
				return(1);
			}
			break;
	  case 'r':
		minR2 = atof(optarg);
		if(minR2 < 0){
			std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a negative minimum R-squared value" << std::endl;
			return(1);
		} else if(minR2 > 1){
			std::cout << Tomahawk::Helpers::timestamp("ERROR")<< "Cannot have minimum R-squared value > 1" << std::endl;
			return(1);
		}
		break;

	  case 'R':
		maxR2 = atof(optarg);
		if(maxR2 < 0){
			std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a negative maximum R-squared value" << std::endl;
		return(1);
		} else if(maxR2 > 1){
			std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a maximum R-squared value > 1" << std::endl;
		return(1);
		}
		break;

	  case 'p':
		  phased = true;
		  forceFunction = true;
		  break;

	  case 'u':
		  phased = false;
		  forceFunction = true;
		  break;

	  case 'P':
		  minP = atof(optarg);
		  if(minP < 0){
			  std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a negative cutoff P-value" << std::endl;
			return(1);
		  } else if(minP > 1){
			  std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a cutoff P-value > 1" << std::endl;
			return(1);
		  }
		  break;
	  case 'a':
		minAlleles = atoi(optarg);
		if(minAlleles < 0){
			std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have negative minimum allele count" << std::endl;
			return(1);
		}
		break;

	  case 'A':
		  maxAlleles = atoi(optarg);
		if(maxAlleles < 0){
			std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have negative maximum allele count" << std::endl;
			return(1);
		}
		break;

	  case 'w':
		  windowBases = atoi(optarg);
		if(windowBases <= 0){
			std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a non-positive window size" << std::endl;
			return(1);
		}
		break;

	  case 'W':
		  windowPosition = atoi(optarg);
		if(windowPosition <= 0){
			std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Cannot have a non-positive window size" << std::endl;
			return(1);
		}
		break;

	  case 's':
		  silent = true;
		  break;

	  default:
		  std::cout << Tomahawk::Helpers::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
		  return(1);
		}
	}

	std::cerr << "here" << std::endl;

	if(input.length() == 0){
		std::cout << Tomahawk::Helpers::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(output.length() == 0){
		std::cout << Tomahawk::Helpers::timestamp("ERROR") << "No output value specified..." << std::endl;
		return(1);
	}

	// Print messages
	programMessage();
	std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling calc..." << std::endl;

	// Parse Totempole index
	std::string index = input + '.' + Tomahawk::Constants::OUTPUT_INDEX_SUFFIX;
	Tomahawk::TotempoleReader totempole;
	if(!totempole.Open(index)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Failed build!" << std::endl;
		return 1;
	}

	// Parse Tomahawk
	Tomahawk::TomahawkReader tomahawk(totempole);
	if(!tomahawk.Open(input)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Failed build!" << std::endl;
		return 1;
	}

	// Parse Tomahawk data
	if(!tomahawk.ValidateHeader()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Failed build!" << std::endl;
		return 1;
	}

	//std::vector<U32> rets = totempole.findOverlaps(Tomahawk::Interval(0, 2221297, 10108169));
	//std::cerr << "Found overlaps: " << rets.size() << std::endl;
	//for(U32 i = 0; i < rets.size(); ++i)
	//	std::cerr << rets[i] << std::endl;

	/*
	std::vector<U32> blocks;
	for(U32 i = 10; i < 12; ++i)
		blocks.push_back(i);

	std::vector< std::pair<U32, U32> > blocks2;
	blocks2.push_back(std::pair<U32,U32>(0, 10));
	*/

	//blocks.push_back(19);
	//blocks.push_back(20);
	//blocks.push_back(21);
	//tomahawk.getBlocks(blocks);
	if(!tomahawk.SetPThreshold(minP))
		return false;

	if(!tomahawk.SetR2Threshold(minR2,maxR2))
		return false;

	tomahawk.SetMinimumAlleles(minAlleles);
	if(!tomahawk.SetThreads(threads))
		return false;

	tomahawk.SetSilent(silent);
	if(forceFunction)
		tomahawk.SetPhased(phased);

	tomahawk.SetOutputType(type);
	if(!tomahawk.OpenWriter(output))
		return false;

	if(!tomahawk.SetChunkDesired(parts))
		return false;

	if(!tomahawk.SetChunkSelected(startPart))
		return false;

	tomahawk.Calculate();
	//tomahawk.Calculate(blocks2);
	//tomahawk.AllVersusAll(blocks);

	return 0;
}
