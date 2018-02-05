#ifndef TOMAHAWKENTRYMETA_H_
#define TOMAHAWKENTRYMETA_H_

namespace Tomahawk{
// Size of meta entry BEFORE run entries
#define TOMAHAWK_ENTRY_META_SIZE	(sizeof(U32) + sizeof(BYTE) + 2*sizeof(float))

/*
 TomahawkEntryMetaBase is used for reinterpreting
 byte streams of TomahawkEntryMeta
 Number of runs can be inferred from the sample
 number and byte length of the stream
 */
struct TomahawkEntryMetaBase{
	typedef TomahawkEntryMetaBase self_type;

public:
	TomahawkEntryMetaBase() :
		missing(0),
		phased(0),
		position(0),
		ref_alt(0),
		MAF(0),
		HWE_P(0)
	{}

	TomahawkEntryMetaBase(const char* const buffer_stream) :
		missing(0),
		phased(0),
		position(0),
		ref_alt(*reinterpret_cast<const BYTE* const>(&buffer_stream[sizeof(U32)])),
		MAF(*reinterpret_cast<const float* const>(&buffer_stream[sizeof(U32)+sizeof(BYTE)])),
		HWE_P(*reinterpret_cast<const float* const>(&buffer_stream[sizeof(U32)+sizeof(BYTE)+sizeof(float)]))
	{
		//std::cerr << "in base ctor: " << std::endl;
		// Overflow fill packed bits
		//U32* t = reinterpret_cast<U32*>(this->missing);
		//*t = *reinterpret_cast<const U32* const>(&buffer_stream[0]);
		memcpy(this, buffer_stream, sizeof(U32));
		//std::cerr << this->missing << ',' << this->phased << ',' << this->position << '\t' << (int)this->ref_alt << '\t' << this->MAF << std::endl;
	}

	~TomahawkEntryMetaBase(){}

	bool isSingleton(void) const{ return(this->MAF == 0); }

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' << (int)entry.ref_alt << '\t' << entry.MAF << '\t' << entry.HWE_P;
		return(out);
	}

public:
	U32 missing: 1, phased: 1, position: 30;
	BYTE ref_alt;
	float MAF;
	float HWE_P;
};

/*
 TomahawkEntryMeta encodes for the basic information
 regaring a variant line such as position, if any genotypes
 are missing and if the data is phased.
 */
template <class T>
struct TomahawkEntryMeta : public TomahawkEntryMetaBase{
	typedef TomahawkEntryMeta     self_type;
	typedef TomahawkEntryMetaBase parent_type;

public:
	TomahawkEntryMeta() :
		runs(0)
	{}

	// Copy from stream
	TomahawkEntryMeta(const char* const buffer_stream) :
		parent_type(buffer_stream),
		runs(*reinterpret_cast<const T* const>(&buffer_stream[TOMAHAWK_ENTRY_META_SIZE]))
	{
		//std::cerr << "runs: " << (int)this->runs << std::endl;
	}

	~TomahawkEntryMeta(){}

	inline bool isValid(void) const{ return(this->runs > 0); }

	friend std::ofstream& operator<<(std::ofstream& os, const self_type& entry){
		const U32 writePos = ( entry.position << 30 ) | ( entry.phased << 1 ) | entry.missing;
		os.write(reinterpret_cast<const char*>(&writePos), sizeof(U32));
		os.write(reinterpret_cast<const char*>(&entry.ref_alt), sizeof(BYTE));
		os.write(reinterpret_cast<const char*>(&entry.runs), sizeof(U32));
		os.write(reinterpret_cast<const char*>(&entry.MAF), sizeof(float));
		os.write(reinterpret_cast<const char*>(&entry.HWE_P), sizeof(float));
		return os;
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' << (int)entry.ref_alt << '\t' << entry.runs << '\t' << entry.MAF << '\t' << entry.HWE_P;
		return(out);
	}

public:
	T runs; // number of runs
};

}

#endif /* TOMAHAWKENTRYMETA_H_ */
