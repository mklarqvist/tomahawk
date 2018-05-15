#ifndef TOMAHAWKENTRYMETA_H_
#define TOMAHAWKENTRYMETA_H_

#include "support/MagicConstants.h"

namespace tomahawk{

#define TOMAHAWK_ENTRY_META_SIZE (sizeof(U32) + sizeof(BYTE) + sizeof(double)*2 + sizeof(U32))

struct MetaEntry{
	typedef MetaEntry self_type;

public:
	MetaEntry() :
		has_missing(0),
		all_phased(0),
		position(0),
		ref_alt(0),
		AF(0),
		HWE_P(0),
		runs(0)
	{}

	MetaEntry(const char* const buffer_stream) :
		position(*reinterpret_cast<const U32* const>(&buffer_stream[0]) >> 2),
		all_phased((*reinterpret_cast<const U32* const>(&buffer_stream[0]) >> 1) & 1),
		has_missing(*reinterpret_cast<const U32* const>(&buffer_stream[0]) & 1),
		ref_alt(*reinterpret_cast<const BYTE* const>(&buffer_stream[sizeof(U32)])),
		AF(*reinterpret_cast<const double* const>(&buffer_stream[sizeof(U32)+sizeof(BYTE)])),
		HWE_P(*reinterpret_cast<const double* const>(&buffer_stream[sizeof(U32)+sizeof(BYTE)+sizeof(double)])),
		runs(*reinterpret_cast<const U32* const>(&buffer_stream[sizeof(U32)+sizeof(BYTE)+sizeof(double)*2]))
	{

	}
	~MetaEntry() = default;

	inline bool isSingleton(void) const{ return(this->AF == 0); }
	inline const char getRefAllele(void) const{ return(constants::REF_ALT_LOOKUP[this->ref_alt >> 4]); }
	inline const char getAltAllele(void) const{ return(constants::REF_ALT_LOOKUP[this->ref_alt & ((1 << 4) - 1)]); }
	inline const char getPhaseVCFCharacter(void) const{ return(this->all_phased == 1 ? '|' : '/'); }

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const MetaEntry& meta_entry){
		U32 pstring = meta_entry.position << 2;
		pstring |= meta_entry.all_phased << 1;
		pstring |= meta_entry.has_missing;
		buffer  += pstring;
		buffer  += meta_entry.ref_alt;
		buffer  += meta_entry.AF;
		buffer  += meta_entry.HWE_P;
		buffer  += meta_entry.runs;
		return(buffer);
	}

public:
	U32 has_missing: 1, // record have >0 missing values
	    all_phased:  1,  // all genotypes are phased
	    position:    30;   // genomic position
	BYTE ref_alt;       // packed reference/alt information
	double AF;
	double HWE_P;
	U32    runs; // number of runs
};

}

#endif /* TOMAHAWKENTRYMETA_H_ */
