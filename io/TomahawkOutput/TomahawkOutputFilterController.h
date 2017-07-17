#ifndef TOMAHAWKOUTPUTFILTERCONTROLLER_H_
#define TOMAHAWKOUTPUTFILTERCONTROLLER_H_

#include <iostream>
#include <iomanip>

#include "../../TypeDefinitions.h"
#include "TomahawkOutputEntry.h"
#include "../../tomahawk/MagicConstants.h"

namespace Tomahawk {

struct TomahawkOutputFilterRegion{
	typedef TomahawkOutputFilterRegion filter_type;
public:
	U32 position;
	U32 contigID;
	filter_type* linked;
};

class TomahawkOutputFilterController {
	typedef TomahawkOutputFilterController self_type;

public:
	TomahawkOutputFilterController();
	~TomahawkOutputFilterController();

	inline void trigger(void){ this->any_filter_user_set = true; }
	inline const bool& isAnySet(void) const{ return(this->any_filter_user_set); }

	bool setFilterTable(const S32& a, const S32& b, const S32& c, const S32& d){
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

	bool setFilterTable(const S32& all){
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

	bool setFilterD(const float& min, const float& max){
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

	bool setFilterDprime(const float& min, const float& max){
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

	bool setFilterRsquared(const float& min, const float& max){
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

	bool setFilterP(const double& min, const double& max){
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

	bool setFilterPmodel(const double& min, const double& max){
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

	bool setFilterChiSquared(const double& min, const double& max){
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

	bool setFilterMGF(const float& minMGF){
		if(minMGF < 0){
			std::cerr << "impossible" << std::endl;
			return false;
		}
		this->minMGF = minMGF;
		this->trigger();
		return true;
	}

	void setFilterInclude(const U16& val){ this->filterValueInclude = val; this->trigger(); }
	void setFilterExclude(const U16& val){ this->filterValueExclude = val; this->trigger(); }

	bool filter(const IO::TomahawkOutputEntry& target) const;
	bool filterAlleleCount(const IO::TomahawkOutputEntry& target) const;

private:
	template <class T, class Y> bool filter(const T& value, const Y& min, const Y& max) const;

private:
	bool any_filter_user_set;
	U64 filteredCount;
	U64 keepCount;

	// Filters
	U32 minP1, minP2, minQ1, minQ2;
	U32 minMinorAlleles, maxMinorAlleles;
	float minD, maxD;
	float minDprime, maxDprime;
	double minR2, maxR2;
	double minP, maxP;
	double minChiSquared, maxChiSquared;
	double minPmodel, maxPmodel;
	float minMGF;

	U16 filterValueInclude;
	U16 filterValueExclude;

	// Todo: filter region tree
};

template <class T, class Y>
inline bool TomahawkOutputFilterController::filter(const T& value, const Y& min, const Y& max) const{
	if(value < min || value > max)
		return false;

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTFILTERCONTROLLER_H_ */
