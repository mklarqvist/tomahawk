#ifndef TOMAHAWKENTRYMETA_H_
#define TOMAHAWKENTRYMETA_H_

namespace Tomahawk{
namespace Support{

// Size of meta entry BEFORE run entries
#define TOMAHAWK_ENTRY_META_SIZE	(sizeof(BYTE) + sizeof(U32) + sizeof(BYTE) + 2*sizeof(float) + sizeof(U32))

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
	TomahawkEntryMetaBase() :
		missing(0),
		phased(0),
		biallelic(0),
		simple(0),
		unused(0),
		position(0),
		ref_alt(0),
		MGF(0),
		HWE_P(0),
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
		BYTE controller = 0;
		controller |= entry.simple << 3;
		controller |= entry.biallelic << 2;
		controller |= entry.phased << 1;
		controller |= entry.missing;
		buffer += controller;
		buffer += entry.position;
		buffer += entry.ref_alt;
		buffer += entry.MGF;
		buffer += entry.HWE_P;
		buffer += entry.virtual_offset;
		return(buffer);
	}

public:
	BYTE missing: 1, phased: 1, biallelic: 1, simple: 1, unused: 4;
	U32 position;
	BYTE ref_alt; // most sites are bi-allelic SNV
	// MGF and HWE is pre-computed as they are used in
	// LD selection
	float MGF;
	float HWE_P;
	// offset to the byte end of this entry in the stream
	U32 virtual_offset;
};

/*
 TomahawkEntryMeta encodes for the basic information
 regarding a variant line such as position, if any genotypes
 are missing and if the all the data is phased.
 */
#pragma pack(1)
template <class T>
struct TomahawkEntryMeta : public TomahawkEntryMetaBase{
	typedef TomahawkEntryMeta self_type;

public:
	TomahawkEntryMeta() : runs(0){}
	~TomahawkEntryMeta(){}

	inline bool isValid(void) const{ return(this->runs > 0); }

	friend std::ofstream& operator<<(std::ofstream& os, const self_type& entry){
		const U32 writePos = ( entry.position << 29 ) | ( entry.biallelic << 2 ) << ( entry.phased << 1 ) | entry.missing;
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
