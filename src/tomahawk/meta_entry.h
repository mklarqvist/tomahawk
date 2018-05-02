#ifndef TOMAHAWKENTRYMETA_H_
#define TOMAHAWKENTRYMETA_H_

namespace Tomahawk{
// Size of meta entry BEFORE run entries
#define TOMAHAWK_ENTRY_META_SIZE (sizeof(U32) + sizeof(BYTE) + 2*sizeof(double))

/*
 TomahawkEntryMetaBase is used for reinterpreting
 byte streams of TomahawkEntryMeta
 Number of runs can be inferred from the sample
 number and byte length of the stream
 */
#pragma pack(push, 1)
struct __attribute__((packed, aligned(1))) MetaEntryBase{
	typedef MetaEntryBase self_type;

public:
	MetaEntryBase() :
		has_missing(0),
		all_phased(0),
		position(0),
		ref_alt(0),
		AF(0),
		HWE_P(0)
	{}

	MetaEntryBase(const char* const buffer_stream){ memcpy(this, buffer_stream, sizeof(self_type)); }
	~MetaEntryBase() = default;

	inline bool isSingleton(void) const{ return(this->AF == 0); }
	inline const char getRefAllele(void) const{ return(Constants::REF_ALT_LOOKUP[this->ref_alt >> 4]); }
	inline const char getAltAllele(void) const{ return(Constants::REF_ALT_LOOKUP[this->ref_alt & ((1 << 4) - 1)]); }
	inline const char getPhaseVCFCharacter(void) const{ return(this->all_phased == 1 ? '|' : '/'); }

	// Overloaded operator for debug use
	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' << (int)entry.ref_alt << '\t' << entry.AF << '\t' << entry.HWE_P;
		return(out);
	}

public:
	U32 has_missing: 1, // record have >0 missing values
	    all_phased:  1,  // all genotypes are phased
	    position:    30;   // genomic position
	BYTE ref_alt;       // packed reference/alt information
	double AF;
	double HWE_P;
};

/*
 TomahawkEntryMeta encodes for the basic information
 regaring a variant line such as position, if any genotypes
 are missing and if the data is phased.
 */
template <class T>
struct __attribute__((packed, aligned(1))) MetaEntry : public MetaEntryBase{
	typedef MetaEntry     self_type;
	typedef MetaEntryBase parent_type;

public:
	MetaEntry() :
		runs(0)
	{}

	// Copy from stream
	MetaEntry(const char* const buffer_stream) :
		parent_type(buffer_stream),
		runs(*reinterpret_cast<const T* const>(&buffer_stream[TOMAHAWK_ENTRY_META_SIZE]))
	{
		//std::cerr << "runs: " << (int)this->runs << std::endl;
	}

	~MetaEntry(){}

	inline bool isValid(void) const{ return(this->runs > 0); }

	friend std::ofstream& operator<<(std::ofstream& os, const self_type& entry){
		const U32 writePos = ( entry.position << 30 ) | ( entry.all_phased << 1 ) | entry.has_missing;
		os.write(reinterpret_cast<const char*>(&writePos),      sizeof(U32));
		os.write(reinterpret_cast<const char*>(&entry.ref_alt), sizeof(BYTE));
		os.write(reinterpret_cast<const char*>(&entry.runs),    sizeof(U32));
		os.write(reinterpret_cast<const char*>(&entry.AF),     sizeof(double));
		os.write(reinterpret_cast<const char*>(&entry.HWE_P),   sizeof(double));
		return os;
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' << (int)entry.ref_alt << '\t' << entry.runs << '\t' << entry.AF << '\t' << entry.HWE_P;
		return(out);
	}

public:
	T runs; // number of runs
};
#pragma pack(pop)

}

#endif /* TOMAHAWKENTRYMETA_H_ */
