#ifndef TOMAHAWK_TOMAHAWKSUPPORT_H_
#define TOMAHAWK_TOMAHAWKSUPPORT_H_

#include "support/MagicConstants.h"

namespace tomahawk{
namespace support{

#pragma pack(push, 1)
template <class T>
struct __attribute__((packed, aligned(1))) GenotypeDiploidRun{
public:
	GenotypeDiploidRun(){}
	GenotypeDiploidRun(const char* const buffer){
		T* t = reinterpret_cast<T*>(this->alleleA);
		*t   = *reinterpret_cast<T*>(buffer);
	}
	~GenotypeDiploidRun(){}

	T alleleA: constants::TOMAHAWK_ALLELE_PACK_WIDTH,
	  alleleB: constants::TOMAHAWK_ALLELE_PACK_WIDTH,
	  runs:    sizeof(T)*8 - constants::TOMAHAWK_SNP_PACK_WIDTH;
};


template <class T>
struct __attribute__((packed, aligned(1))) GenotypeDiploidRunPacked{
public:
	GenotypeDiploidRunPacked(){}
	GenotypeDiploidRunPacked(const char* const buffer){
		T* t = reinterpret_cast<T*>(this->alleles);
		*t   = *reinterpret_cast<T*>(buffer);
	}
	~GenotypeDiploidRunPacked(){}

	T alleles: constants::TOMAHAWK_SNP_PACK_WIDTH,
	  runs:    sizeof(T)*8 - constants::TOMAHAWK_SNP_PACK_WIDTH;
};

#pragma pack(pop)

} // end support

namespace Constants{

// Dummy variables used for Twk > VCF output
const std::string PASS("PASS");
const std::string GT("GT");
const std::string QUAL("100");

}

}



#endif /* TOMAHAWK_TOMAHAWKSUPPORT_H_ */
