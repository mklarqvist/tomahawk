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
		explicit __meta_controller(void) :
				missing(0),
				phased(0),
				biallelic(0),
				simple(0),
				hasComplex(0),
				rle(0),
				unused(0)
		{}
		~__meta_controller(){}

		BYTE missing: 1,  // any missing
		     phased: 1,   // all phased
			 biallelic: 1,// is biallelic
			 simple: 1,   // is simple SNV
			 hasComplex: 1,// has complex meta
			 rle: 1,      // uses RLE compression
			 unused: 2;
	} controller_byte;

public:
	TomahawkEntryMetaBase() :
		position(0),
		ref_alt(0),
		MGF(0),
		HWE_P(0),
		AF(0),
		virtual_offset_complex(0),
		virtual_offset(0)
	{}
	~TomahawkEntryMetaBase(){}

	inline const bool isSingleton(void) const{ return(this->MGF == 0); }
	inline const bool isSimpleSNV(void) const{ return(this->controller.biallelic == true && this->controller.simple == true); }
	inline const bool isRLE(void) const{ return(this->controller.rle); }

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' << (int)entry.controller.biallelic << ',' << (int)entry.controller.simple << '\t' << (int)entry.ref_alt << '\t' << entry.MGF << '\t' << entry.HWE_P << '\t' << entry.virtual_offset;
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
	// MGF, HWE, AF is pre-computed as they are used in
	// LD heuristics
	float MGF;
	float HWE_P;
	float AF;
	// Hot-cold split structure. pointer to cold data
	// since a pointer cannot be read from a byte
	// stream as its memory location changes
	// we have to provide the pointer as an ABSOLUTE
	// virtual stream offset relative to the complex
	// start position into the complex byte stream
	U32 virtual_offset_complex;
	// offset to the byte end of this entry in the stream
	// (either RLE or simple stream depending on context)
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
struct TomahawkEntryMeta : public TomahawkEntryMetaBase{
	typedef TomahawkEntryMeta self_type;

public:
	TomahawkEntryMeta() : n_runs(0){}
	~TomahawkEntryMeta(){}

	inline const bool isValid(void) const{ return(this->n_runs > 0); }
	inline const T& getRuns(void) const{ return(this->n_runs); }

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' << (int)entry.controller.biallelic << ',' << (int)entry.controller.simple << '\t' << (int)entry.ref_alt << '\t' << entry.n_runs << '\t' << entry.MGF << '\t' << entry.HWE_P << '\t' << entry.n_runs;
		return(out);
	}

public:
	T n_runs; // number of runs
};

}
}

#endif /* TOMAHAWKENTRYMETA_H_ */
