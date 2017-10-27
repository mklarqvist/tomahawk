#ifndef TOMAHAWK_TOMAHAWKGTENTRIES_H_
#define TOMAHAWK_TOMAHAWKGTENTRIES_H_

namespace Tomahawk{
namespace Support{

// We CANNOT place phasing template parameter
// since if we set it to 0 then we have a
// 0 width bit field which is illegal
// To solve this we introduce the TomahawkRunNoPhase
// data structure below
#pragma pack(1)
template <class T, BYTE missing = 1>
struct TomahawkRun{
private:
	typedef TomahawkRun self_type;

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
	  alleleA: (1 + missing),
	  alleleB: (1 + missing),
	  runs:    sizeof(T)*8 - (2 * (1 + missing) + 1);
};

#pragma pack(1)
template <class T, BYTE missing = 1>
struct TomahawkRunNoPhase{
private:
	typedef TomahawkRunNoPhase self_type;

public:
	TomahawkRunNoPhase();  // Disallowed ctor
	~TomahawkRunNoPhase(); // Disallowed dtor

	T alleleA: (1 + missing),
	  alleleB: (1 + missing),
	  runs:    sizeof(T)*8 - (2 * (1 + missing));
};

#pragma pack(1)
template <class T, BYTE missing = 1>
struct TomahawkRunPacked{
private:
	typedef TomahawkRunPacked self_type;

public:
	TomahawkRunPacked();	// Disallowed ctor
	~TomahawkRunPacked(); // Disallowed dtor

	T phasing: 1,
	  alleles: 2*(1 + missing),
	  runs:    sizeof(T)*8 - (2 * (1 + missing) + 1);
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
