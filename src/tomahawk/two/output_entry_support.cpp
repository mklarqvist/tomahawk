#include "../two/output_entry_support.h"

namespace Tomahawk {
namespace Support {

OutputEntrySupport::OutputEntrySupport() :
	controller(0),
	R2(0),
	D(0),
	Dprime(0),
	Dmax(0),
	P(0),
	chiSqModel(0),
	chiSqFisher(0),
	totalAlleleCounts(0)
	{
	}

OutputEntrySupport::~OutputEntrySupport(){ }


void OutputEntrySupport::operator=(const OutputEntrySupport& other){
	this->R2 = other.R2;
	this->D = other.D;
	this->Dprime = other.Dprime;
	this->P = other.P;
	this->totalAlleleCounts = other.totalAlleleCounts;
	memcpy(&this->alleleCounts[0], &other.alleleCounts[0], sizeof(float)*171);
	memcpy(&this->haplotypeCounts[0], &other.haplotypeCounts[0], sizeof(float)*4);
}

void OutputEntrySupport::printUnphasedCounts(void) const{
	// Prints 3x3 Punnett square for unphased data
	std::cerr << this->alleleCounts[0] << '\t' << this->alleleCounts[1] + this->alleleCounts[4] << '\t' << this->alleleCounts[5] << '\t'
			  << this->alleleCounts[16] + this->alleleCounts[64] << '\t' << this->alleleCounts[17] + this->alleleCounts[20] + this->alleleCounts[65] + this->alleleCounts[68] << '\t' << this->alleleCounts[21] + this->alleleCounts[69] << '\t'
			  << this->alleleCounts[80] << '\t' << this->alleleCounts[81]+this->alleleCounts[84] << '\t' << this->alleleCounts[85] << std::endl;
}

void OutputEntrySupport::printPhasedCounts(void) const{
	std::cerr << this->alleleCounts[0] << '\t' << this->alleleCounts[1] << '\t' << this->alleleCounts[4] << '\t' << this->alleleCounts[5] << std::endl;
}

}
} /* namespace Tomahawk */
