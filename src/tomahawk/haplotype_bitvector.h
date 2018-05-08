#ifndef TOMAHAWK_HAPLOTYPE_BITVECTOR_H_
#define TOMAHAWK_HAPLOTYPE_BITVECTOR_H_

namespace Tomahawk{
namespace Base{

struct HaplotypeBitVector{
public:
	HaplotypeBitVector() :
		n_bytes(0),
		l_list(0),
		indices(nullptr),
		entries(nullptr)
	{

	}

	HaplotypeBitVector(const U64 n_entries) :
		n_bytes(ceil((double)n_entries/64)),
		l_list(0),
		indices(nullptr),
		entries(new U64[n_bytes])
	{
		memset(this->entries, 0, this->n_bytes);
	}
	~HaplotypeBitVector(){ delete [] this->entries; delete [] this->indices; }

	inline void reset(void){ memset(this->entries, 0, this->n_bytes); }
	inline const bool operator[](const U32& position) const{ return(this->entries[position/64] & (1 << (position % 64))); }
	inline const bool get(const U32& position) const{ return(this->entries[position/64] & (1 << (position % 64)));}
	inline void set(const U32& position){ this->entries[position/64] |= (1 << (position % 64)); }
	inline void set(const U32& position, const bool val){ this->entries[position/64] |= (val << (position % 64)); }

public:
	U32 n_bytes;
	U32 l_list; // number of odd items
	U32* indices;
	U64* entries;
};

}
}

#endif /* TOMAHAWK_HAPLOTYPE_BITVECTOR_H_ */
