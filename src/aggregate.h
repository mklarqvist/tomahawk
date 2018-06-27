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
#ifndef AGGREGATE_H_
#define AGGREGATE_H_


#include "tomahawk/tomahawk_output_reader.h"
#include "utility.h"

void aggregate_usage(void){
	programMessage();
	std::cerr <<
	"About:  Aggregate terrabytes of TWO data into rasterized objects of size [x,y] for\n"
	"        plotting. Data has to be sorted and pre-sliced to the correct interval-of-interest\n"
	"        prior to running."
	"Usage:  " << tomahawk::constants::PROGRAM_NAME << " sort [options] <in.two>\n\n"
	"Options:\n"
	"  -i FILE       input Tomahawk output file (TWO; required)\n"
	"  -x INT        X-axis scene dimension in pixels (required)\n"
	"  -y INT        Y-axis scene dimension in pixels (required)\n"
	//"  -X INT-INT    X-axis data range (from, to)\n"
	//"  -Y INT-INT    Y-axis data range (from, to)\n"
	"  -d STRING     Data field to aggregate: May be one of {R, R2, D, DPrime, p1, p2, q1, q2, P}(required)\n"
	"  -f STRING     Aggregation function: May be one of {count, min, max, mean, sd, sum, sum-squared}(required)\n"
	//"  -T STRING     Transformation function: May be one of {linear, squared, sqrt}(required)\n"
	"  -s            Hide all program messages [null]\n\n";
}

int aggregate(int argc, char** argv){
	if(argc < 3){
		aggregate_usage();
		return(1);
	}

	static struct option long_options[] = {
		{"input",           required_argument, 0, 'i' },
		{"scene-x",         required_argument, 0, 'x' },
		{"scene-y",         required_argument, 0, 'y' },
		{"target-field",    required_argument, 0, 'd' },
		{"reduction-field", required_argument, 0, 'f' },
		{"silent",          no_argument,       0, 's' },
		{0,0,0,0}
	};

	// Parameter defaults
	std::string input, temp;
	tomahawk::support::aggregation_parameters parameters;

	S32 scene_x_pixels     = -1;
	S32 scene_y_pixels     = -1;
	S32 aggregation_target = -1;
	S32 reduction_target   = -1;

	int c = 0;
	int long_index = 0;
	while ((c = getopt_long(argc, argv, "i:x:y:X:Y:d:f:s", long_options, &long_index)) != -1){
		switch (c){
		case ':':   /* missing option argument */
			fprintf(stderr, "%s: option `-%c' requires an argument\n",
					argv[0], optopt);
			break;

		case 'i':
			input = std::string(optarg);
			break;

		case 'x':
			parameters.scene_x_pixels = atoi(optarg);
			scene_x_pixels = atoi(optarg);
			break;

		case 'y':
			parameters.scene_y_pixels = atoi(optarg);
			scene_y_pixels = atoi(optarg);
			break;

		case 'd':
			temp = std::string(optarg);
			if(strncasecmp("R", temp.data(), 1) == 0 && temp.size() == 1){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_R;
			} else if(strncasecmp("R2", temp.data(), 2) == 0 && temp.size() == 2){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_R2;
			} else if(strncasecmp("D", temp.data(), 1) == 0 && temp.size() == 1){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_D;
			} else if(strncasecmp("Dprime", temp.data(), 6) == 0 && temp.size() == 6){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_DPrime;
			} else if(strncasecmp("P", temp.data(), 1) == 0 && temp.size() == 1){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_P_VALUE;
			} else if(strncasecmp("P1", temp.data(), 2) == 0 && temp.size() == 2){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_P1;
			} else if(strncasecmp("P2", temp.data(), 2) == 0 && temp.size() == 2){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_P2;
			} else if(strncasecmp("Q1", temp.data(), 2) == 0 && temp.size() == 2){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_Q1;
			} else if(strncasecmp("Q2", temp.data(), 2) == 0 && temp.size() == 2){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_Q2;
			} else if(strncasecmp("logP", temp.data(), 4) == 0 && temp.size() == 4){
				parameters.aggregation_target = tomahawk::support::TWK_AGGREGATE_LOG_P_VALUE;
			} else {
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Unknown aggregation target: " << temp << "..." << std::endl;
				return(1);
			}
			aggregation_target = parameters.aggregation_target;
			break;

		case 'f':
			// count, min, max, mean, sd, sum, sum-squared
			temp = std::string(optarg);
			if(strncasecmp("count", temp.data(), 5) == 0 && temp.size() == 5){
				parameters.reduction_target = tomahawk::support::TWK_AGGREGATE_REDUCE_COUNT;
			} else if(strncasecmp("min", temp.data(), 3) == 0 && temp.size() == 3){
				parameters.reduction_target = tomahawk::support::TWK_AGGREGATE_REDUCE_MIN;
			} else if(strncasecmp("max", temp.data(), 3) == 0 && temp.size() == 3){
				parameters.reduction_target = tomahawk::support::TWK_AGGREGATE_REDUCE_MAX;
			} else if(strncasecmp("sd", temp.data(), 2) == 0 && temp.size() == 2){
				parameters.reduction_target = tomahawk::support::TWK_AGGREGATE_REDUCE_SD;
			} else if(strncasecmp("sum", temp.data(), 3) == 0 && temp.size() == 3){
				parameters.reduction_target = tomahawk::support::TWK_AGGREGATE_REDUCE_SUM;
			} else if(strncasecmp("sum-squared", temp.data(), 11) == 0 && temp.size() == 11){
				parameters.reduction_target = tomahawk::support::TWK_AGGREGATE_REDUCE_SUM_SQUARED;
			} else if(strncasecmp("mean", temp.data(), 4) == 0 && temp.size() == 4){
				parameters.reduction_target = tomahawk::support::TWK_AGGREGATE_REDUCE_MEAN;
			} else {
				std::cerr << tomahawk::helpers::timestamp("ERROR") << "Unknown reduction function: " << temp << "..." << std::endl;
				return(1);
			}
			reduction_target = parameters.reduction_target;
			break;

		case '?':
		default:
			fprintf(stderr, "%s: option `-%c' is invalid: ignored\n", argv[0], optopt);
			break;
		}
	}

	if(input.length() == 0){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "No input file specified..." << std::endl;
		return(1);
	}

	if(scene_x_pixels == -1){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "Requires scene-X dimension..." << std::endl;
		return(1);
	}

	if(scene_y_pixels == -1){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "Requires scene-Y dimension..." << std::endl;
		return(1);
	}

	if(reduction_target == -1){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "Requires reduction target..." << std::endl;
		return(1);
	}

	if(aggregation_target == -1){
		std::cerr << tomahawk::helpers::timestamp("ERROR") << "Requires aggregation target..." << std::endl;
		return(1);
	}

	if(!SILENT){
		programMessage();
		std::cerr << tomahawk::helpers::timestamp("LOG") << "Calling aggregation..." << std::endl;
	}

	tomahawk::TomahawkOutputReader reader;
	if(!reader.open(input))  return 1;
	if(!reader.aggregate(parameters)) return 1;


	return 0;
}



#endif /* AGGREGATE_H_ */