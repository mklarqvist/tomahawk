#ifndef TOMAHAWKENTRYMETA_H_
#define TOMAHAWKENTRYMETA_H_

namespace Tomahawk{
// Size of meta entry BEFORE run entries
#define TOMAHAWHK_ENTRY_META_SIZE	(sizeof(U32) + sizeof(BYTE) + sizeof(double) + sizeof(double))

#pragma pack(1)
template <class T>
struct TomahawkEntryMeta{
public:
	TomahawkEntryMeta() : position(0), missing(0), phased(0), ref_alt(0), runs(0), MAF(0), HWE_P(0){}
	~TomahawkEntryMeta(){}

	inline bool isValid(void) const{ return(this->runs > 0); }

	friend std::ofstream& operator<<(std::ofstream& os, const TomahawkEntryMeta& entry){
		const U32 writePos = ( entry.position << 30 ) | ( entry.phased << 1 ) | entry.missing;
		os << writePos << entry.ref_alt << entry.runs << entry.MAF << entry.HWE_P;
		return os;
	}

	friend std::ostream& operator<<(std::ostream& out, const TomahawkEntryMeta& entry){
		out << (entry.position >> 2) << '\t' << (int)entry.ref_alt << '\t' << entry.runs << '\t' << entry.MAF << '\t' << entry.HWE_P;
		return(out);
	}

public:
	// Todo: change doubles to float
	// memory aligned to 16 byte boundaries
	const U32 missing: 1, phased: 1, position: 30;
	const BYTE ref_alt;
	const double MAF;
	const double HWE_P;
	T runs;
};

}

#endif /* TOMAHAWKENTRYMETA_H_ */
