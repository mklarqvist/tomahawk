#ifndef TOMAHAWKOUTPUTFILTERCONTROLLER_H_
#define TOMAHAWKOUTPUTFILTERCONTROLLER_H_

#include <iostream>
#include <iomanip>
#include <vector>

#include "../../support/TypeDefinitions.h"
#include "TomahawkOutputEntry.h"
#include "../../support/MagicConstants.h"

namespace Tomahawk {

class TomahawkOutputFilterController {
	typedef TomahawkOutputFilterController self_type;
	typedef IO::TomahawkOutputEntry entry_type;
	typedef bool (self_type::*filterFunction)(const entry_type& entry) const;

public:
	TomahawkOutputFilterController();
	~TomahawkOutputFilterController();

	inline const bool& isAnySet(void) const{ return(this->any_filter_user_set); }

	bool setFilterTable(const S32& a, const S32& b, const S32& c, const S32& d);
	bool setFilterTable(const S32& all);
	bool setFilterD(const float& min, const float& max);
	bool setFilterDprime(const float& min, const float& max);
	bool setFilterRsquared(const float& min, const float& max);
	bool setFilterP(const double& min, const double& max);
	bool setFilterPmodel(const double& min, const double& max);
	bool setFilterChiSquared(const double& min, const double& max);
	void setFilterInclude(const U16& val){ this->filterValueInclude = val; this->trigger(); }
	void setFilterExclude(const U16& val){ this->filterValueExclude = val; this->trigger(); }
	bool setFilterJointHF(const int64_t& min, const int64_t& max);

	bool filter(const entry_type& target) const;

private:
	inline void trigger(void){ this->any_filter_user_set = true; }

	template <class T, class Y> bool filter(const T& value, const Y& min, const Y& max) const;

	bool filterJointHF(const entry_type& target) const;
	bool filterHF(const entry_type& target) const;
	bool filterD(const entry_type& type) const;
	bool filterDprime(const entry_type& type) const;
	bool filterRsquared(const entry_type& type) const;
	bool filterP(const entry_type& type) const;
	bool filterPmodel(const entry_type& type) const;
	bool filterFLAG(const entry_type& type) const;

private:
	bool any_filter_user_set;

	// Filters
	U64 minP1, minP2, minQ1, minQ2;
	U64 maxP1, maxP2, maxQ1, maxQ2;
	U64 minMHF, maxMHF;
	float minD, maxD;
	float minDprime, maxDprime;
	double minR2, maxR2;
	double minP, maxP;
	double minChiSquared, maxChiSquared;
	double minPmodel, maxPmodel;
	U16 filterValueInclude;
	U16 filterValueExclude;

	std::vector<filterFunction> filter_functions; // push filter functions to array and loop over
};

template <class T, class Y>
inline bool TomahawkOutputFilterController::filter(const T& value, const Y& min, const Y& max) const{
	if(value < min || value > max)
		return false;

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTFILTERCONTROLLER_H_ */
