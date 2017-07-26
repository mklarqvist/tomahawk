#include <limits>
#include "TomahawkOutputFilterController.h"

namespace Tomahawk {

TomahawkOutputFilterController::TomahawkOutputFilterController() :
	any_filter_user_set(false),
	filteredCount(0),
	keepCount(0),
	minP1(std::numeric_limits<U64>::max()),
	minP2(std::numeric_limits<U64>::max()),
	minQ1(std::numeric_limits<U64>::max()),
	minQ2(std::numeric_limits<U64>::max()),
	minMHF(0),
	maxMHF(std::numeric_limits<U64>::max()),
	minD(-100), maxD(100),
	minDprime(-100), maxDprime(100),
	minR2(0), maxR2(100),
	minP(0), maxP(100),
	minChiSquared(0), maxChiSquared(std::numeric_limits<double>::max()),
	minPmodel(0), maxPmodel(1),
	filterValueInclude(0),
	filterValueExclude(0)
{}

TomahawkOutputFilterController::~TomahawkOutputFilterController(){}

bool TomahawkOutputFilterController::filter(const IO::TomahawkOutputEntry& target) const{
	if(((target.FLAGS & this->filterValueInclude) != this->filterValueInclude) || ((target.FLAGS & this->filterValueExclude) != 0)){
		//std::cerr << "failed bits" << std::endl;
		return false;
	}
	if(!this->filter(target.R2, this->minR2, this->maxR2)) return false;
	if(!this->filter(target.D, this->minD, this->maxD)) return false;
	if(!this->filter(target.P, this->minP, this->maxP)) return false;
	if(!this->filter(target.Dprime, this->minDprime, this->maxDprime)) return false;
	if(!this->filter(target.chiSqModel, this->minChiSquared, this->maxChiSquared)) return false;
	if(!this->filterAlleleCount(target)) return false;

	return true;
}

bool TomahawkOutputFilterController::filterAlleleCount(const IO::TomahawkOutputEntry& target) const{
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

} /* namespace Tomahawk */
