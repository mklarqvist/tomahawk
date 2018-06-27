#ifndef TOMAHAWKOUTPUTFILTERCONTROLLER_H_
#define TOMAHAWKOUTPUTFILTERCONTROLLER_H_

#include <iostream>
#include <iomanip>
#include <vector>

#include "../../support/magic_constants.h"
#include "support/type_definitions.h"
#include "output_entry.h"

namespace tomahawk {

class OutputFilter {
	typedef OutputFilter self_type;
	typedef io::OutputEntry entry_type;
	typedef bool (self_type::*filterFunction)(const entry_type& entry) const;

public:
	OutputFilter();
	~OutputFilter();

	inline const bool& isAnySet(void) const{ return(this->any_filter_user_set); }

	bool setFilterTable(const double& a, const double& b, const double& c, const double& d);
	bool setFilterTable(const double& all);
	inline void setFilterUpperTriangular(const bool set){ this->upper_triangular_only = set; }
	inline void setFilterLowerTriangular(const bool set){ this->lower_triangular_only = set; }
	bool setFilterD(const float& min, const float& max);
	bool setFilterDprime(const float& min, const float& max);
	bool setFilterR(const float& min, const float& max);
	bool setFilterRsquared(const float& min, const float& max);
	bool setFilterP(const double& min, const double& max);
	bool setFilterPmodel(const double& min, const double& max);
	bool setFilterChiSquaredTable(const double& min, const double& max);
	bool setFilterChiSquaredModel(const double& min, const double& max);

	void setFilterInclude(const U16& val){ this->FLAGInclude = val; this->trigger(); }
	void setFilterExclude(const U16& val){ this->FLAGExclude = val; this->trigger(); }
	bool setFilterMHF(const double& min, const double& max);

	bool filter(const entry_type& target) const;

	std::string getInterpretedString(void) const;
	inline void trigger(void){ this->any_filter_user_set = true; }

private:
	template <class T, class Y> bool filter(const T& value, const Y& min, const Y& max) const;

	bool filterJointHF(const entry_type& target) const;
	bool filterHF(const entry_type& target) const;

public:
	bool any_filter_user_set;
	bool upper_triangular_only;
	bool lower_triangular_only;

	// Filters
	double minP1, minP2, minQ1, minQ2;
	double maxP1, maxP2, maxQ1, maxQ2;
	double minMHF, maxMHF;
	float  minD, maxD;
	float  minDprime, maxDprime;
	double minR, maxR;
	double minR2, maxR2;
	double minP, maxP;
	double minChiSquaredTable, maxChiSquaredTable;
	double minChiSquaredModel, maxChiSquaredModel;
	U16    FLAGInclude;
	U16    FLAGExclude;

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
