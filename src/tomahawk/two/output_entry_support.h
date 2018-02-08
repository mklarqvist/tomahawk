#ifndef TOMAHAWK_TOMAHAWKOUTPUTLD_H_
#define TOMAHAWK_TOMAHAWKOUTPUTLD_H_

#include "../../io/BasicBuffer.h"

namespace Tomahawk {
namespace Support{

struct OutputEntrySupport{
private:
	typedef OutputEntrySupport self_type;

public:
	OutputEntrySupport();
	~OutputEntrySupport();

	inline void resetPhased(void){
		this->controller      = 0;
		this->alleleCounts[0] = 0;
		this->alleleCounts[1] = 0;
		this->alleleCounts[4] = 0;
		this->alleleCounts[5] = 0;
		this->chiSqModel      = 0;
		// All other values can legally overflow
		// They are not used
	}

	inline void resetUnphased(void){
		this->controller       = 0;
		this->alleleCounts[0]  = 0;
		this->alleleCounts[1]  = 0;
		this->alleleCounts[4]  = 0;
		this->alleleCounts[5]  = 0;
		this->alleleCounts[16] = 0;
		this->alleleCounts[17] = 0;
		this->alleleCounts[20] = 0;
		this->alleleCounts[21] = 0;
		this->alleleCounts[64] = 0;
		this->alleleCounts[65] = 0;
		this->alleleCounts[68] = 0;
		this->alleleCounts[69] = 0;
		this->alleleCounts[80] = 0;
		this->alleleCounts[81] = 0;
		this->alleleCounts[84] = 0;
		this->alleleCounts[85] = 0;
		// All other values can legally overflow
		// They are not used
	}

	inline float& operator[](const U32& p){ return(this->alleleCounts[p]); }
	inline const float& operator[](const U32& p) const{ return(this->alleleCounts[p]); }
	void operator=(const OutputEntrySupport& other);
	void printUnphasedCounts(void) const;
	void printPhasedCounts(void) const;

	friend IO::BasicBuffer& operator<<(IO::BasicBuffer& os, const OutputEntrySupport& entry){
		// Notice that CONTROLLER is written separately
		os += entry.alleleCounts[0];
		os += entry.alleleCounts[1];
		os += entry.alleleCounts[4];
		os += entry.alleleCounts[5];
		os += entry.D;
		os += entry.Dprime;
		os += entry.R2;
		os += entry.P;
		os += entry.chiSqFisher;
		os += entry.chiSqModel;
		return os;
	}

	friend std::ostream& operator<<(std::ostream& os, const OutputEntrySupport& entry){
		os << entry.alleleCounts[0] << '\t' << entry.alleleCounts[1] << '\t' << entry.alleleCounts[4] << '\t' << entry.alleleCounts[5] << '\t'
		   << entry.D << '\t' << entry.Dprime << '\t' << entry.R2 << '\t' << entry.P << '\t' << entry.chiSqFisher;
		return os;
	}

	inline void setPhased(void)          { this->controller |= 1;   }
	inline void setHasMissingValues(void){ this->controller |= 2;   }
	inline void setIncomplete(void)      { this->controller |= 4;   }
	inline void setMultipleRoots(void)   { this->controller |= 8;   }
	inline void setSameContig(void)      { this->controller |= 16;  }
	inline void setLongRange(void)       { this->controller |= 32;  }
	inline void setFailedHWEA(void)      { this->controller |= 64;  }
	inline void setFailedHWEB(void)      { this->controller |= 128; }
	inline void setLowMAFA(void)         { this->controller |= 256; }
	inline void setLowMAFB(void)         { this->controller |= 512; }

	inline const float countAlternatives(void) const{
		// Find largest
		const float* max = &this->alleleCounts[0];
		if(this->alleleCounts[1] > *max) max = &this->alleleCounts[1];
		if(this->alleleCounts[4] > *max) max = &this->alleleCounts[4];
		if(this->alleleCounts[5] > *max) max = &this->alleleCounts[5];

		// Count number of non-major allele combinations
		float max2 = 0;
		if(&this->alleleCounts[0] != max) max2 += this->alleleCounts[0];
		if(&this->alleleCounts[1] != max) max2 += this->alleleCounts[1];
		if(&this->alleleCounts[4] != max) max2 += this->alleleCounts[4];
		if(&this->alleleCounts[5] != max) max2 += this->alleleCounts[5];

		return(max2);
	}

public:
	U16    controller;        // FLAG byte
	float  R2;                // R squared
	float  D;                 // D
	float  Dprime;            // D'
	float  Dmax;              // Dmax
	double P;                 // Fisher or Chi-Squared P value for 2x2 contingency table
	double chiSqModel;        // Chi-Squared critical value for 3x3 contingency table
	double chiSqFisher;       // Chi-Squared critical value for 2x2 contingency table
	float  totalAlleleCounts; // Total number of alleles

	// Counters
	float  alleleCounts[171];
	float  haplotypeCounts[4];
};

} /* namespace Support */
} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKOUTPUTLD_H_ */
