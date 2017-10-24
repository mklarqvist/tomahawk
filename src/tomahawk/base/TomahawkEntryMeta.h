#ifndef TOMAHAWKENTRYMETA_H_
#define TOMAHAWKENTRYMETA_H_

#include "../../io/BasicBuffer.h"

namespace Tomahawk{
namespace Support{

// Size of meta entry BEFORE run entries
#define TOMAHAWK_ENTRY_META_SIZE	(sizeof(BYTE) + sizeof(U32) + sizeof(BYTE) + 2*sizeof(float) + 2*sizeof(U32))

/*
 TomahawkEntryMetaBase is used for reinterpreting
 byte streams of TomahawkEntryMeta
 Number of runs can be inferred from the sample
 number and byte length of the stream
 */
#pragma pack(1)
struct TomahawkEntryMetaBase{
	typedef TomahawkEntryMetaBase self_type;
	typedef IO::BasicBuffer buffer_type;

public:
	typedef struct __meta_controller{
		explicit __meta_controller(void) : missing(0), phased(0), biallelic(0), simple(0), hasComplex(0), unused(0){}
		~__meta_controller(){}

		BYTE missing: 1, phased: 1, biallelic: 1, simple: 1, hasComplex: 1, unused: 3;
	} controller_byte;

public:
	TomahawkEntryMetaBase() :
		position(0),
		ref_alt(0),
		MGF(0),
		HWE_P(0),
		virtual_offset_complex(0),
		virtual_offset(0)
	{}
	~TomahawkEntryMetaBase(){}

	inline bool isSingleton(void) const{ return(this->MGF == 0); }

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' << (int)entry.ref_alt << '\t' << entry.MGF << '\t' << entry.HWE_P << '\t' << entry.virtual_offset;
		return(out);
	}

	// Overload operator+= for basic buffer
	friend buffer_type& operator+=(buffer_type& buffer, const self_type& entry){
		buffer += *reinterpret_cast<const BYTE* const>(&entry.controller);
		buffer += entry.position;
		buffer += entry.ref_alt;
		buffer += entry.MGF;
		buffer += entry.HWE_P;
		buffer += entry.virtual_offset_complex;
		buffer += entry.virtual_offset;
		return(buffer);
	}

public:
	controller_byte controller;
	U32 position;
	// most sites are bi-allelic and simple SNVs
	// sites that are not bi-allelic and not simple
	// will be encoded in the complex meta section
	BYTE ref_alt;
	// MGF and HWE is pre-computed as they are used in
	// LD heuristics
	float MGF;
	float HWE_P;
	// Hot-cold split structure. pointer to cold data
	// since a pointer cannot be read from a byte
	// stream as its memory location changes
	// we have to provide the pointer as an ABSOLUTE
	// virtual stream offset relative to the complex
	// start position into the complex byte stream
	U32 virtual_offset_complex;
	// offset to the byte end of this entry in the stream
	// this allows fast iteration when switching between
	// the two compression approaches
	U32 virtual_offset;
};

/*
 TomahawkEntryMetaRLE encodes for the basic information
 regarding a variant line such as position, if any genotypes
 are missing and if the all the data is phased.
 */
#pragma pack(1)
template <class T>
struct TomahawkEntryMetaRLE : public TomahawkEntryMetaBase{
	typedef TomahawkEntryMetaRLE self_type;

public:
	TomahawkEntryMetaRLE() : runs(0){}
	~TomahawkEntryMetaRLE(){}

	inline bool isValid(void) const{ return(this->runs > 0); }

	friend std::ofstream& operator<<(std::ofstream& os, const self_type& entry){
		const U32 writePos = ( entry.position << 29 ) | ( entry.controller.biallelic << 2 ) << ( entry.controller.phased << 1 ) | entry.controller.missing;
		os.write(reinterpret_cast<const char*>(&writePos),      sizeof(U32));
		os.write(reinterpret_cast<const char*>(&entry.ref_alt), sizeof(BYTE));
		os.write(reinterpret_cast<const char*>(&entry.runs),    sizeof(U32));
		os.write(reinterpret_cast<const char*>(&entry.MGF),     sizeof(float));
		os.write(reinterpret_cast<const char*>(&entry.HWE_P),   sizeof(float));
		return os;
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' << (int)entry.ref_alt << '\t' << entry.runs << '\t' << entry.MGF << '\t' << entry.HWE_P;
		return(out);
	}

public:
	T runs; // number of runs
};

}
}

#endif /* TOMAHAWKENTRYMETA_H_ */
