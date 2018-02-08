#ifndef TOMAHAWKOUTPUTENTRY_H_
#define TOMAHAWKOUTPUTENTRY_H_

#include "../../io/BasicBuffer.h"
#include "../../totempole/TotempoleContig.h"

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
	typedef OutputEntry                    self_type;
	typedef Totempole::TotempoleContigBase contig_type;

	// if interpreted directly from buffer stream
	OutputEntry()  = default;
	~OutputEntry() = default;

	// Copy data from stream
	OutputEntry(const char* const data_buffer){
		memcpy(this, data_buffer, sizeof(self_type));
	}

	// Copy data from stream
	OutputEntry(const IO::BasicBuffer& data_buffer){
		memcpy(this, data_buffer.data, sizeof(self_type));
	}

	OutputEntry(const self_type* const other){
		memcpy(this, other, sizeof(self_type));
	}

	// Comparator function
	// Called from sort helper only
	bool operator<(const self_type& other) const{
		if (this->AcontigID < other.AcontigID) return true;
		if (other.AcontigID < this->AcontigID) return false;

		if (this->Aposition < other.Aposition) return true;
		if (other.Aposition < this->Aposition) return false;

		if (this->BcontigID < other.BcontigID) return true;
		if (other.BcontigID < this->BcontigID) return false;

		if (this->Bposition < other.Bposition) return true;
		if (other.Bposition < this->Bposition) return false;

		return false;
	}

	// Comparator function: inverse of lesser comparator
	bool operator>(const self_type& other){ return(!((*this) < other)); }

	friend std::ostream& operator<<(std::ostream& os, const self_type& entry){
		os << std::setprecision(8) << (int)entry.FLAGS << '\t' << entry.AcontigID << '\t' << entry.Aposition << '\t' << entry.BcontigID << '\t' << entry.Bposition
			<< '\t' << entry.p1 << '\t' << entry.p2 << '\t' << entry.q1 << '\t' << entry.q2 << '\t' << entry.D << '\t' << entry.Dprime
			<< '\t' << entry.R2 << '\t' << entry.P << '\t' << entry.chiSqFisher << '\t' << entry.chiSqModel;

		return(os);
	}

	std::ostream& write(std::ostream& os, const contig_type* const contigs) const{
		os << std::setprecision(8) << (int)this->FLAGS << '\t' << contigs[this->AcontigID].name << '\t' << this->Aposition << '\t' << contigs[this->BcontigID].name << '\t' << this->Bposition
			<< '\t' << this->p1 << '\t' << this->p2 << '\t' << this->q1 << '\t' << this->q2 << '\t' << this->D << '\t' << this->Dprime
			<< '\t' << this->R2 << '\t' << this->P << '\t' << this->chiSqFisher << '\t' << this->chiSqModel << '\n';

		return(os);
	}

	friend IO::BasicBuffer& operator<<(IO::BasicBuffer& b, const self_type& entry){
		b.Add(reinterpret_cast<const char*>(&entry), sizeof(self_type));
		return(b);
	}

	// Swaps cA,pA with cB,pB
	// used in sorting for indices
	void swapDirection(void){
		U32 Ac = this->AcontigID;
		U32 Bc = this->BcontigID;
		this->AcontigID = Bc;
		this->BcontigID = Ac;
		U32& A = *reinterpret_cast<U32*>(((char*)this + sizeof(U16) +   sizeof(U32)));
		U32& B = *reinterpret_cast<U32*>(((char*)this + sizeof(U16) + 3*sizeof(U32)));
		std::swap(A,B);
	}

	U16    FLAGS;
	U32    AcontigID;
	U32    Amissing: 1, Aphased: 1, Aposition: 30;
	U32    BcontigID;
	U32    Bmissing: 1, Bphased: 1, Bposition: 30;
	float  p1, p2, q1, q2;
	float  D, Dprime;
	float  R2;
	double P;
	double chiSqFisher;
	double chiSqModel;
};


// Sort reinterpreted casts of data
// Workaround is based on data being reinterpreted as
// entries from byte streams. Therefore regular sorting
// is illegal. Instead a BYTE array literals stored as
// a hard copy and reinterpreted as an entry in the
// overloaded operator<
struct __attribute__((packed, aligned(1))) OutputEntrySort{
	typedef OutputEntrySort self_type;
	typedef OutputEntry     parent_type;

	OutputEntrySort(){}
	OutputEntrySort(const self_type& other){
		memcpy(this->data, other.data, sizeof(parent_type));
	}
	OutputEntrySort(self_type&& other) noexcept{ std::swap(this->data, other.data); }
	OutputEntrySort& operator=(const self_type& other){
		self_type tmp(other);	// re-use copy-constructor
		*this = std::move(tmp); // re-use move-assignment
		return *this;
	}
	OutputEntrySort& operator=(self_type&& other) noexcept{
		std::swap(this->data, other.data);
		return *this;
	}
	~OutputEntrySort(){ }

	bool operator<(const self_type& other) const{
		const parent_type& self_parent  = *reinterpret_cast<const parent_type* const>(&this->data[0]);
		const parent_type& other_parent = *reinterpret_cast<const parent_type* const>(&other.data[0]);
		return(self_parent < other_parent);
	}

	BYTE data[sizeof(parent_type)]; // reinterpret me as entry
};
#pragma pack(pop)

// comparator functions for output entry
namespace Support{

static inline bool TomahawkOutputEntryCompFuncConst(const OutputEntry& self, const OutputEntry& other){
	if (self.AcontigID < other.AcontigID) return true;
	if (other.AcontigID < self.AcontigID) return false;

	if (self.Aposition < other.Aposition) return true;
	if (other.Aposition < self.Aposition) return false;

	if (self.BcontigID < other.BcontigID) return true;
	if (other.BcontigID < self.BcontigID) return false;

	if (self.Bposition < other.Bposition) return true;
	if (other.Bposition < self.Bposition) return false;

	return false;
}

}

}
}

#endif /* TOMAHAWKOUTPUTENTRY_H_ */
