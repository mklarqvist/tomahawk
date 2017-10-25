#ifndef TOMAHAWK_TOMAHAWKGTENTRIES_H_
#define TOMAHAWK_TOMAHAWKGTENTRIES_H_

namespace Tomahawk{
namespace Support{

#pragma pack(1)
template <class T>
struct TomahawkRun{
private:
	typedef TomahawkRun<T> self_type;

public:
	TomahawkRun();	// Disallowed ctor
	~TomahawkRun(); // Disallowed dtor

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << (U32)entry.alleleA
			<< (entry.phasing ? '|' : '/')
			<< (U32)entry.alleleB;
		return(out);
	}

	T phasing: 1,
	  alleleA: Constants::TOMAHAWK_ALLELE_PACK_WIDTH,
	  alleleB: Constants::TOMAHAWK_ALLELE_PACK_WIDTH,
	  runs:    sizeof(T)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH - 1;
};

#pragma pack(1)
template <class T>
struct TomahawkRunPacked{
public:
	TomahawkRunPacked();	// Disallowed ctor
	~TomahawkRunPacked();	// Disallowed dtor

	T phasing: 1,
	  alleles: Constants::TOMAHAWK_SNP_PACK_WIDTH,
	  runs:    sizeof(T)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH - 1;
};

#pragma pack(1)
template <class T>
struct TomahawkRunSimple{
private:
	typedef TomahawkRunSimple<T> self_type;

public:
	TomahawkRunSimple();  // Disallowed ctor
	~TomahawkRunSimple(); // Disallowed dtor

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << (entry.alleleA == 0 ? '.' : (U16)entry.alleleA - 1)
			<< (entry.phasing ? '|' : '/')
			<< (entry.alleleB == 0 ? '.' : (U16)entry.alleleB - 1);
		return(out);
	}

	T phasing: 1,
	  alleleA: (sizeof(T)*8)/2 - 1,
	  alleleB: (sizeof(T)*8)/2 - 1,
	  unused:  1;
};

}
}

#endif /* TOMAHAWK_TOMAHAWKGTENTRIES_H_ */
