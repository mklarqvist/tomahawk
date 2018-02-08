#include <limits>
#include <string>
#include "output_filter.h"

namespace Tomahawk {

OutputFilter::OutputFilter() :
	any_filter_user_set(false),
	minP1(0),
	minP2(0),
	minQ1(0),
	minQ2(0),
	maxP1(std::numeric_limits<double>::max()),
	maxP2(std::numeric_limits<double>::max()),
	maxQ1(std::numeric_limits<double>::max()),
	maxQ2(std::numeric_limits<double>::max()),
	minMHF(0),
	maxMHF(std::numeric_limits<double>::max()),
	minD(-100), maxD(100),
	minDprime(-100), maxDprime(100),
	minR2(0), maxR2(100),
	minP(0), maxP(100),
	minChiSquared(0), maxChiSquared(std::numeric_limits<double>::max()),
	minPmodel(0), maxPmodel(std::numeric_limits<double>::max()),
	filterValueInclude(0),
	filterValueExclude(0)
{}

OutputFilter::~OutputFilter(){}

std::string OutputFilter::getInterpretedString(void) const{
	if(this->any_filter_user_set){
		return(std::string(
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
				"minP=" + std::to_string(this->minP) + " " +
				"maxP=" + std::to_string(this->maxP) + " " +
				"minChiSquared=" + std::to_string(this->minChiSquared) + " " +
				"maxChiSquared=" + (this->maxChiSquared == std::numeric_limits<double>::max() ? "DOUBLE_MAX" : std::to_string(this->maxChiSquared)) + " " +
				"minPmodel=" + std::to_string(this->minPmodel) + " " +
				"maxPmodel=" + std::to_string(this->maxPmodel) + " " +
				"filterValueInclude=" + std::to_string(this->filterValueInclude) + " " +
				"filterValueExclude=" + std::to_string(this->filterValueExclude)
		));
	} else {
		return(std::string("no_filter"));
	}
}

bool OutputFilter::filter(const entry_type& target) const{
	if(((target.FLAGS & this->filterValueInclude) != this->filterValueInclude) || ((target.FLAGS & this->filterValueExclude) != 0))
		return false;

	if(!this->filter(target.R2, this->minR2, this->maxR2)) return false;
	if(!this->filter(target.P, this->minP, this->maxP)) return false;
	if(!this->filter(target.D, this->minD, this->maxD)) return false;
	if(!this->filter(target.Dprime, this->minDprime, this->maxDprime)) return false;
	if(!this->filter(target.chiSqModel, this->minChiSquared, this->maxChiSquared)) return false;
	if(!this->filter(target.chiSqFisher, this->minPmodel, this->maxPmodel)) return false;
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
	const float* max = &target.p1;
	if(target.p2 > *max) max = &target.p2;
	if(target.q1 > *max) max = &target.q1;
	if(target.q2 > *max) max = &target.q2;

	// sum of cells excluding largest
	float total = 0;
	if(&target.p1 != max) total += target.p1;
	if(&target.p2 != max) total += target.p2;
	if(&target.q1 != max) total += target.q1;
	if(&target.q2 != max) total += target.q2;

	return(total > this->minMHF);
}

bool OutputFilter::setFilterTable(const double& a, const double& b, const double& c, const double& d){
	if(a < 0 || b < 0 || c < 0 || d < 0){
		std::cerr << "cannot have negative filter values" << std::endl;
		return false;
	}

	if(a == 0 && b == 0 && c == 0 && d == 0){
		std::cerr << "cannot filter with all cells set to 0" << std::endl;
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
		std::cerr << "cannot have negative filter values" << std::endl;
		return false;
	}
	if(all == 0){
		std::cerr << "cannot filter with all cells set to 0" << std::endl;
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

	this->minD = min - Constants::ALLOWED_ROUNDING_ERROR;
	this->maxD = max + Constants::ALLOWED_ROUNDING_ERROR;
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

	this->minDprime = min - Constants::ALLOWED_ROUNDING_ERROR;
	this->maxDprime = max + Constants::ALLOWED_ROUNDING_ERROR;
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

	this->minR2 = min - Constants::ALLOWED_ROUNDING_ERROR;
	this->maxR2 = max + Constants::ALLOWED_ROUNDING_ERROR;
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

bool OutputFilter::setFilterPmodel(const double& min, const double& max){
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

	this->minPmodel = min;
	this->maxPmodel = max;
	this->trigger();
	return true;
}

bool OutputFilter::setFilterChiSquared(const double& min, const double& max){
	if(max < min){
		std::cerr << "max < min" << std::endl;
		return false;
	}
	if(min < 0 || max < 0){
		std::cerr << "has negative" << std::endl;
		return false;
	}

	this->minChiSquared = min;
	this->maxChiSquared = max;
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
