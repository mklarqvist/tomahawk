#ifndef TOMAHAWK_HAPLOTYPE_BITVECTOR_H_
#define TOMAHAWK_HAPLOTYPE_BITVECTOR_H_

namespace tomahawk{
namespace base{

struct HaplotypeBitVector{
public:
	HaplotypeBitVector() :
		n_bytes(0),
		n_actual(0),
		l_list(0),
		indices(nullptr),
		sampled_indicies(nullptr),
		entries(nullptr)
	{

	}

	HaplotypeBitVector(const U64 n_entries) :
		n_bytes(ceil((double)n_entries/8)),
		n_actual(0),
		l_list(0),
		indices(nullptr),
		sampled_indicies(nullptr),
		entries(new BYTE[n_bytes])
	{
		memset(this->entries, 0, this->n_bytes);
	}
	~HaplotypeBitVector(){
		delete [] this->entries;
		delete [] this->indices;
		delete [] this->sampled_indicies;
	}

	inline void reset(void){ memset(this->entries, 0, this->n_bytes); }
	inline const bool operator[](const U32& position) const{ return(this->entries[position/8] & (1 << (position % 8))); }
	inline const bool get(const U32& position) const{ return(this->entries[position/8] & (1 << (position % 8)));}
	inline void set(const U32& position){ this->entries[position/8] |= (1 << (position % 8)); }
	inline void set(const U32& position, const bool val){ this->entries[position/8] |= (val << (position % 8)); }

public:
	U32 n_bytes;
	U32 n_actual;
	U32 l_list; // number of odd items
	U32* indices;
	U32* sampled_indicies;
	BYTE* entries;
};

}
}

#endif /* TOMAHAWK_HAPLOTYPE_BITVECTOR_H_ */
