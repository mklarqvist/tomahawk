#ifndef TOMAHAWKOUTPUTENTRY_H_
#define TOMAHAWKOUTPUTENTRY_H_

#include "../../io/BasicBuffer.h"
#include "../../index/index_contig.h"

namespace Tomahawk{
namespace IO{

/**<
 * Primary data structure for Tomahawk-generated LD output.
 * This higher-order primitive can be interpreted directly
 * from a packed buffer (unaligned memory access) or
 * explicitly by invoking a valid constructor
 */
#pragma pack(push, 1)
struct __attribute__((packed, aligned(1))) OutputEntry{
public:
	typedef OutputEntry                    self_type;
	typedef Totempole::HeaderContig        contig_type;
	typedef IO::BasicBuffer                buffer_type;

public:
	// if interpreted directly from buffer stream
	OutputEntry() :
		FLAGS(0),
		AcontigID(0),
		Amissing(false),
		Aphased(false),
		Aposition(0),
		BcontigID(0),
		Bmissing(false),
		Bphased(false),
		Bposition(0),
		p1(0), p2(0), q1(0), q2(0),
		D(0), Dprime(0),
		R(0), R2(0),
		P(0),
		chiSqFisher(0),
		chiSqModel(0)
	{

	}

	~OutputEntry() = default;
	// Copy data from stream
	OutputEntry(const char* const data_buffer);
	// Copy data from stream
	OutputEntry(const buffer_type& data_buffer);
	OutputEntry(const self_type* const other);

	// Comparator function
	// Called from sort helper only
	inline static bool sortDescending(const self_type& a, const self_type& b){ return(a < b); }
	inline static bool sortAscending(const self_type& a, const self_type& b){ return(a > b); }
	bool operator< (const self_type& other) const;
	bool operator<=(const self_type& other) const;
	bool operator==(const self_type& other) const;

	// Comparator function: inverse of lesser comparator
	inline bool operator> (const self_type& other) const{ return(!((*this) <  other)); }
	inline bool operator>=(const self_type& other) const{ return(!((*this) <= other)); }

	// Swaps cA,pA with cB,pB
	// used in sorting for indices
	void swapDirection(void);

	friend std::ostream& operator<<(std::ostream& os, const self_type& entry){
		os << std::setprecision(8) << (int)entry.FLAGS << '\t' << entry.AcontigID << '\t' << entry.Aposition << '\t' << entry.BcontigID << '\t' << entry.Bposition
			<< '\t' << entry.p1 << '\t' << entry.p2 << '\t' << entry.q1 << '\t' << entry.q2 << '\t' << entry.D << '\t' << entry.Dprime
			<< '\t' << entry.R << '\t' << entry.R2 << '\t' << entry.P << '\t' << entry.chiSqFisher << '\t' << entry.chiSqModel;

		return(os);
	}

	std::ostream& write(std::ostream& os, const contig_type* const contigs) const{
		os << std::setprecision(8) << (int)this->FLAGS << '\t' << contigs[this->AcontigID].name << '\t' << this->Aposition << '\t' << contigs[this->BcontigID].name << '\t' << this->Bposition
			<< '\t' << this->p1 << '\t' << this->p2 << '\t' << this->q1 << '\t' << this->q2 << '\t' << this->D << '\t' << this->Dprime
			<< '\t' << this->R << '\t' << this->R2 << '\t' << this->P << '\t' << this->chiSqFisher << '\t' << this->chiSqModel << '\n';

		return(os);
	}

	// Write to buffer
	friend buffer_type& operator<<(buffer_type& b, const self_type& entry){
		b.Add(reinterpret_cast<const char*>(&entry), sizeof(self_type));
		return(b);
	}

public:
	U16    FLAGS;
	U32    AcontigID;
	U32    Amissing: 1, Aphased: 1, Aposition: 30;
	U32    BcontigID;
	U32    Bmissing: 1, Bphased: 1, Bposition: 30;
	float  p1, p2, q1, q2;
	float  D, Dprime; // D and D'
	float  R, R2;     // Correlation coefficient
	double P;         // P-value
	double chiSqFisher;
	double chiSqModel;
};
#pragma pack(pop)

}
}

#endif /* TOMAHAWKOUTPUTENTRY_H_ */
