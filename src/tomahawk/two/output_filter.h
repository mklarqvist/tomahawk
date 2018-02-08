#ifndef TOMAHAWKOUTPUTFILTERCONTROLLER_H_
#define TOMAHAWKOUTPUTFILTERCONTROLLER_H_

#include <iostream>
#include <iomanip>
#include <vector>

#include "../../support/type_definitions.h"
#include "../../support/MagicConstants.h"
#include "../two/output_entry.h"

namespace Tomahawk {

class OutputFilter {
	typedef OutputFilter self_type;
	typedef IO::OutputEntry entry_type;
	typedef bool (self_type::*filterFunction)(const entry_type& entry) const;

public:
	OutputFilter();
	~OutputFilter();

	inline const bool& isAnySet(void) const{ return(this->any_filter_user_set); }

	bool setFilterTable(const double& a, const double& b, const double& c, const double& d);
	bool setFilterTable(const double& all);
	bool setFilterD(const float& min, const float& max);
	bool setFilterDprime(const float& min, const float& max);
	bool setFilterRsquared(const float& min, const float& max);
	bool setFilterP(const double& min, const double& max);
	bool setFilterPmodel(const double& min, const double& max);
	bool setFilterChiSquared(const double& min, const double& max);
	void setFilterInclude(const U16& val){ this->filterValueInclude = val; this->trigger(); }
	void setFilterExclude(const U16& val){ this->filterValueExclude = val; this->trigger(); }
	bool setFilterMHF(const double& min, const double& max);

	bool filter(const entry_type& target) const;

	std::string getInterpretedString(void) const;
	inline void trigger(void){ this->any_filter_user_set = true; }

private:
	template <class T, class Y> bool filter(const T& value, const Y& min, const Y& max) const;

	bool filterJointHF(const entry_type& target) const;
	bool filterHF(const entry_type& target) const;
	bool filterD(const entry_type& type) const;
	bool filterDprime(const entry_type& type) const;
	bool filterRsquared(const entry_type& type) const;
	bool filterP(const entry_type& type) const;
	bool filterPmodel(const entry_type& type) const;
	bool filterFLAG(const entry_type& type) const;

public:
	bool any_filter_user_set;

	// Filters
	double minP1, minP2, minQ1, minQ2;
	double maxP1, maxP2, maxQ1, maxQ2;
	double minMHF, maxMHF;
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
inline bool OutputFilter::filter(const T& value, const Y& min, const Y& max) const{
	if(value < min || value > max)
		return false;

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTFILTERCONTROLLER_H_ */
