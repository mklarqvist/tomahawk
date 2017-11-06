#ifndef ALGORITHM_COMPRESSION_RADIXSORTGT_H_
#define ALGORITHM_COMPRESSION_RADIXSORTGT_H_

#include "../../io/bcf/BCFReader.h"

namespace Tomahawk {
namespace Algorithm {

/*
 * This class performs a radix sort on a
 * block of variant lines given they are
 * biallelic diploid.
 */
class RadixSortGT {
	typedef RadixSortGT self_type;
	typedef BCF::BCFReader bcf_reader_type;
	typedef BCF::BCFEntry  bcf_entry_type;

public:
	RadixSortGT(const U64 n_samples);
	~RadixSortGT();

	// Reset does NOT need to cast after each
	// iteration as values are overwritten
	// each cycle
	inline void reset(void){
		this->position = 0;
		memset(this->GT_array, 0, sizeof(BYTE)*n_samples);
		memset(&p_i, 0, sizeof(U32)*9);

		for(U32 i = 0; i < this->n_samples; ++i)
			this->ppa[i] = i;
	}

	// Construct given a reader with a block
	// of BCF entries loaded in it
	bool build(const bcf_reader_type& reader);
	bool update(const bcf_entry_type& entry);

	inline const U64& getSamples(void) const{ return(this->n_samples); }
	inline const U32& size(void) const{ return(this->position); }
	inline const U32& operator[](const U32& p) const{return(this->ppa[p]); }
	inline const U32* getPPA(void) const{return(this->ppa); }

	// Debug functions
	void outputGT(const bcf_reader_type& reader);
	bool assesRLECost(const bcf_reader_type& reader);
	U64 assessRLEPPA(const bcf_entry_type& reader);
	U64 assessRLE(const bcf_entry_type& reader);

private:
	U64 n_samples; // total number of entries in file
	U32 position;  // number of entries parsed
	U32 p_i[9];    // number of entries in bin i
	U32* ppa;      // position prefix array
	BYTE* GT_array;// packed genotype array
	U32** bins;    // bin i

	// costs
	U64 cost_ppa_conventional;
	U64 cost_ppa_best;
	U64 cost_ppa_byte;
	U64 cost_ppa_u16;
	U64 cost_ppa_u32;
	U64 cost_ppa_u64;
	U64 cost_rle_conventional;
	U64 cost_rle_best;
	U64 cost_rle_byte;
	U64 cost_rle_u16;
	U64 cost_rle_u32;
	U64 cost_rle_u64;
};

} /* namespace Algorithm */
} /* namespace Tomahawk */

#endif /* ALGORITHM_COMPRESSION_RADIXSORTGT_H_ */
