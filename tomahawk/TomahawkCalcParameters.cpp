#ifndef TOMAHAWK_TOMAHAWKCALCPARAMETERS_CPP_
#define TOMAHAWK_TOMAHAWKCALCPARAMETERS_CPP_

#include "TomahawkCalcParameters.h"

namespace Tomahawk{

bool TomahawkCalcParameters::Validate(void){
	if(n_threads < 0){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Invalid number of threads..." << std::endl;
		return false;
	}

	if(n_chunks < 0){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Invalid number of partitions..." << std::endl;
		return false;
	}

	if(chunk_selected < 0 || chunk_selected > n_chunks){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Invalid selected partition..." << std::endl;
		return false;
	}

	if(R2_min < 0 || R2_min > 1){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Invalid minimum R-squared cutoff..." << std::endl;
		return false;
	}

	if(R2_max < 0 || R2_max > 1){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Invalid maximum R-squared cutoff..." << std::endl;
		return false;
	}

	if(R2_min > R2_max){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Minimum R-squared value > maximum R-squared value..." << std::endl;
		return false;
	}

	if(P_threshold < 0 || P_threshold > 1){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Invalid P-value cutoff..." << std::endl;
		return false;
	}

	if(minimum_alleles < 0){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Invalid minimum number of alleles..." << std::endl;
		return false;
	}

	if(maximum_alleles < 0){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Invalid maximum number of alleles..." << std::endl;
		return false;
	}

	if(minimum_alleles > maximum_alleles){
		std::cerr << Helpers::timestamp("ERROR", "CALC") << "Minimum number of alleles > maximum number of alleles..." << std::endl;
		return false;
	}

	this->R2_min -= Constants::ALLOWED_ROUNDING_ERROR;
	this->R2_max += Constants::ALLOWED_ROUNDING_ERROR;

	return true;
}

}



#endif /* TOMAHAWK_TOMAHAWKCALCPARAMETERS_CPP_ */
