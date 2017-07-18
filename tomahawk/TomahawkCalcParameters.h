#ifndef TOMAHAWK_TOMAHAWKCALCPARAMETERS_H_
#define TOMAHAWK_TOMAHAWKCALCPARAMETERS_H_

#include "../io/BasicWriters.h"

namespace Tomahawk{

// Constant
const double STANDARD_MIN_R2_CUTOFF = 0.3;

struct TomahawkCalcParameters{
	enum force_method {none, phasedFunction, unphasedFunction};

	TomahawkCalcParameters() :
		R2_min(STANDARD_MIN_R2_CUTOFF - Constants::ALLOWED_ROUNDING_ERROR),
		R2_max(1 + Constants::ALLOWED_ROUNDING_ERROR),
		P_threshold(1e-6),
		minimum_alleles(5),
		compression_type(IO::GenericWriterInterace::compression::binary),
		force(force_method::none)
	{}

	TomahawkCalcParameters(const double R2min, const double R2max, const double Pmin, IO::GenericWriterInterace::compression type, force_method force = force_method::none):
		R2_min(R2min),
		R2_max(R2max),
		P_threshold(Pmin),
		minimum_alleles(5),
		compression_type(type),
		force(force)
	{}

	~TomahawkCalcParameters(){}

	friend std::ostream& operator<<(std::ostream& os, const TomahawkCalcParameters& p){
		os << Helpers::timestamp("CALC", "PARAMETERS") << "R-squared (" << p.R2_min << '-' << p.R2_max << "), P < " << p.P_threshold << ", non-refs > " << p.minimum_alleles;
		return(os);
	}

	double R2_min;
	double R2_max;
	double P_threshold;
	U64 minimum_alleles;
	IO::GenericWriterInterace::compression compression_type;
	force_method force;
};

}



#endif /* TOMAHAWK_TOMAHAWKCALCPARAMETERS_H_ */
