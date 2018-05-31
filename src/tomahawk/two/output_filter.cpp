#include <limits>
#include <string>
#include "output_filter.h"

namespace tomahawk {

OutputFilter::OutputFilter() :
	any_filter_user_set(false),
	upper_triangular_only(false),
	lower_triangular_only(false),
	minP1(0),
	minP2(0),
	minQ1(0),
	minQ2(0),
	maxP1(std::numeric_limits<double>::max()),
	maxP2(std::numeric_limits<double>::max()),
	maxQ1(std::numeric_limits<double>::max()),
	maxQ2(std::numeric_limits<double>::max()),
	minMHF(-1),
	maxMHF(std::numeric_limits<double>::max()),
	minD(-100), maxD(100),
	minDprime(-100), maxDprime(100),
	minR(0), maxR(100),
	minR2(0), maxR2(100),
	minP(0), maxP(100),
	minChiSquaredTable(0), maxChiSquaredTable(std::numeric_limits<double>::max()),
	minChiSquaredModel(0), maxChiSquaredModel(std::numeric_limits<double>::max()),
	FLAGInclude(0),
	FLAGExclude(0)
{}

OutputFilter::~OutputFilter(){}

std::string OutputFilter::getInterpretedString(void) const{
	if(this->any_filter_user_set){
		return(std::string(
			"anyFilterSet=" + (this->any_filter_user_set ? std::string("TRUE ") : std::string("FALSE ")) +
			"upperTriangular=" + (this->upper_triangular_only ? std::string("TRUE ") : std::string("FALSE ")) +
			"lowerTriangular=" + (this->upper_triangular_only ? std::string("TRUE ") : std::string("FALSE ")) +
			"minP1=" + std::to_string(this->minP1) + " " +
			"minP2=" + std::to_string(this->minP2) + " " +
			"minQ1=" + std::to_string(this->minQ1) + " " +
			"minQ2=" + std::to_string(this->minQ2) + " " +
			"maxP1=" + (this->maxP1 == std::numeric_limits<double>::max() ? "DOUBLE_MAX" : std::to_string(this->maxP1)) + " " +
			"maxP2=" + (this->maxP2 == std::numeric_limits<double>::max() ? "DOUBLE_MAX" : std::to_string(this->maxP2)) + " " +
			"maxQ1=" + (this->maxQ1 == std::numeric_limits<double>::max() ? "DOUBLE_MAX" : std::to_string(this->maxQ1)) + " " +
			"maxQ2=" + (this->maxQ2 == std::numeric_limits<double>::max() ? "DOUBLE_MAX" : std::to_string(this->maxQ2)) + " " +
			"minMHF=" + std::to_string(this->minMHF) + " " +
			"maxMHF=" + (this->maxMHF == std::numeric_limits<double>::max() ? "DOUBLE_MAX" : std::to_string(this->maxMHF)) + " " +
			"minD=" + std::to_string(this->minD) + " " +
			"maxD=" + std::to_string(this->maxD) + " " +
			"minDprime=" + std::to_string(this->minDprime) + " " +
			"maxDprime=" + std::to_string(this->maxDprime) + " " +
			"minR2=" + std::to_string(this->minR2) + " " +
			"maxR2=" + std::to_string(this->maxR2) + " " +
			"minR=" + std::to_string(this->minR) + " " +
			"maxR=" + std::to_string(this->maxR) + " " +
			"minP=" + std::to_string(this->minP) + " " +
			"maxP=" + std::to_string(this->maxP) + " " +
			"minChiSquaredTable=" + std::to_string(this->minChiSquaredTable) + " " +
			"maxChiSquaredTable=" + (this->maxChiSquaredTable == std::numeric_limits<double>::max() ? "DOUBLE_MAX" : std::to_string(this->maxChiSquaredTable)) + " " +
			"minChiSquaredModel=" + std::to_string(this->minChiSquaredModel) + " " +
			"maxChiSquaredModel=" + (this->maxChiSquaredModel == std::numeric_limits<double>::max() ? "DOUBLE_MAX" : std::to_string(this->maxChiSquaredModel)) + " " +
			"FLAGInclude=" + std::to_string(this->FLAGInclude) + " " +
			"FLAGExclude=" + std::to_string(this->FLAGExclude)
		));
	} else {
		return(std::string("anyFilterSet=FALSE"));
	}
}

bool OutputFilter::filter(const entry_type& target) const{
	if(this->upper_triangular_only){
		if(target.AcontigID == target.BcontigID && target.Aposition > target.Bposition){
			return false;
		}
	}

	if(this->lower_triangular_only){
		if(target.AcontigID == target.BcontigID && target.Aposition < target.Bposition){
			return false;
		}
	}

	if(((target.FLAGS & this->FLAGInclude) != this->FLAGInclude) || ((target.FLAGS & this->FLAGExclude) != 0))
		return false;

	if(!this->filter(target.R, this->minR, this->maxR)) return false;
	if(!this->filter(target.R2, this->minR2, this->maxR2)) return false;
	if(!this->filter(target.P, this->minP, this->maxP)) return false;
	if(!this->filter(target.D, this->minD, this->maxD)) return false;
	if(!this->filter(target.Dprime, this->minDprime, this->maxDprime)) return false;
	if(!this->filter(target.chiSqModel, this->minChiSquaredModel, this->maxChiSquaredModel)) return false;
	if(!this->filter(target.chiSqFisher, this->minChiSquaredTable, this->maxChiSquaredTable)) return false;
	if(!this->filterJointHF(target)) return false;
	if(!this->filter(target.p1, this->minP1, this->maxP1)) return false;
	if(!this->filter(target.p2, this->minP2, this->maxP2)) return false;
	if(!this->filter(target.q1, this->minQ1, this->maxQ1)) return false;
	if(!this->filter(target.q2, this->minQ2, this->maxQ2)) return false;

	return true;
}

bool OutputFilter::filterHF(const entry_type& target) const{
	return(target.p1 >= this->minP1 || target.p2 >= this->minP2 || target.q1 >= this->minQ1 || target.q2 >= this->minQ2);
}

bool OutputFilter::filterJointHF(const entry_type& target) const{
	// find largest
	double max = target.p1;
	BYTE targetField = 0;
	if(target.p2 > max){ max = target.p2; targetField = 1; }
	if(target.q1 > max){ max = target.q1; targetField = 2; }
	if(target.q2 > max){ max = target.q2; targetField = 3; }

	// sum of cells excluding largest
	double total = 0;
	if(targetField != 0) total += target.p1;
	if(targetField != 1) total += target.p2;
	if(targetField != 2) total += target.q1;
	if(targetField != 3) total += target.q2;

	return(total > this->minMHF);
}

bool OutputFilter::setFilterTable(const double& a, const double& b, const double& c, const double& d){
	if(a < 0 || b < 0 || c < 0 || d < 0){
		std::cerr << "Cannot have negative filter values" << std::endl;
		return false;
	}

	if(a == 0 && b == 0 && c == 0 && d == 0){
		std::cerr << "Cannot filter with all cells set to 0" << std::endl;
		return false;
	}

	this->minP1 = a;
	this->minP2 = b;
	this->minQ1 = c;
	this->minQ2 = d;
	this->trigger();
	return true;
}

bool OutputFilter::setFilterTable(const double& all){
	if(all < 0){
		std::cerr << "Cannot have negative filter values" << std::endl;
		return false;
	}
	if(all == 0){
		std::cerr << "Cannot filter with all cells set to 0" << std::endl;
		return false;
	}
	this->minP1 = all;
	this->minP2 = all;
	this->minQ1 = all;
	this->minQ2 = all;
	this->trigger();
	return true;
}

bool OutputFilter::setFilterD(const float& min, const float& max){
	if(max < min){
		std::cerr << "max < min" << std::endl;
		return false;
	}
	if(min < 0 || max < 0){
		std::cerr << "has negative" << std::endl;
		return false;
	}
	if(min > 1 || max > 1){
		std::cerr << "value > 1" << std::endl;
		return false;
	}

	this->minD = min - constants::ALLOWED_ROUNDING_ERROR;
	this->maxD = max + constants::ALLOWED_ROUNDING_ERROR;
	this->trigger();
	return true;
}

bool OutputFilter::setFilterDprime(const float& min, const float& max){
	if(max < min){
		std::cerr << "max < min" << std::endl;
		return false;
	}
	//if(min < 0 || max < 0){
	//	std::cerr << "has negative" << std::endl;
	//	return false;
	//}
	if(min > 1 || max > 1){
		std::cerr << "value > 1" << std::endl;
		return false;
	}

	this->minDprime = min - constants::ALLOWED_ROUNDING_ERROR;
	this->maxDprime = max + constants::ALLOWED_ROUNDING_ERROR;
	this->trigger();
	return true;
}

bool OutputFilter::setFilterR(const float& min, const float& max){
	if(max < min){
		std::cerr << "max < min" << std::endl;
		return false;
	}
	if(min < 0 || max < 0){
		std::cerr << "has negative" << std::endl;
		return false;
	}
	if(min > 1 || max > 1){
		std::cerr << "value > 1" << std::endl;
		return false;
	}

	this->minR = min - constants::ALLOWED_ROUNDING_ERROR;
	this->maxR = max + constants::ALLOWED_ROUNDING_ERROR;
	this->trigger();
	return true;
}


bool OutputFilter::setFilterRsquared(const float& min, const float& max){
	if(max < min){
		std::cerr << "max < min" << std::endl;
		return false;
	}
	if(min < 0 || max < 0){
		std::cerr << "has negative" << std::endl;
		return false;
	}
	if(min > 1 || max > 1){
		std::cerr << "value > 1" << std::endl;
		return false;
	}

	this->minR2 = min - constants::ALLOWED_ROUNDING_ERROR;
	this->maxR2 = max + constants::ALLOWED_ROUNDING_ERROR;
	this->trigger();
	return true;
}

bool OutputFilter::setFilterP(const double& min, const double& max){
	if(max < min){
		std::cerr << "max < min" << std::endl;
		return false;
	}
	if(min < 0 || max < 0){
		std::cerr << "has negative" << std::endl;
		return false;
	}
	if(min > 1 || max > 1){
		std::cerr << "value > 1" << std::endl;
		return false;
	}

	this->minP = min;
	this->maxP = max;
	this->trigger();
	return true;
}

bool OutputFilter::setFilterChiSquaredTable(const double& min, const double& max){
	if(max < min){
		std::cerr << "max < min" << std::endl;
		return false;
	}
	if(min < 0 || max < 0){
		std::cerr << "has negative" << std::endl;
		return false;
	}

	this->minChiSquaredTable = min;
	this->maxChiSquaredTable = max;
	this->trigger();
	return true;
}

bool OutputFilter::setFilterChiSquaredModel(const double& min, const double& max){
	if(max < min){
		std::cerr << "max < min" << std::endl;
		return false;
	}
	if(min < 0 || max < 0){
		std::cerr << "has negative" << std::endl;
		return false;
	}

	this->minChiSquaredModel = min;
	this->maxChiSquaredModel = max;
	this->trigger();
	return true;
}

bool OutputFilter::setFilterMHF(const double& min, const double& max){
	if(min < 0){
		std::cerr << "min < 0" << std::endl;
		return false;
	}

	if(max < min){
		std::cerr << "max < min" << std::endl;
		return false;
	}

	this->minMHF = min;
	this->maxMHF = max;
	this->trigger();
	return true;
}

} /* namespace Tomahawk */
