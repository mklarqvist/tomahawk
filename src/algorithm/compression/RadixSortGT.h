#ifndef ALGORITHM_COMPRESSION_RADIXSORTGT_H_
#define ALGORITHM_COMPRESSION_RADIXSORTGT_H_

namespace Tomahawk {
namespace Algorithm {

class RadixSortGT {
	typedef RadixSortGT self_type;
	typedef BCF::BCFReader bcf_reader_type;
	typedef BCF::BCFEntry  bcf_entry_type;

public:
	RadixSortGT();
	~RadixSortGT();

	void reset(void);

	bool build(bcf_reader_type& reader){
		if(reader.size() == 0)
			return false;

		// Cycle over BCF entries
		for(U32 i = 0; i < reader.size(); ++i){
			// Has to be biallelic
			// otherwise skip
			if(!reader[i].isBiallelic())
				continue;

			if(!this->update(reader[i]))
				continue;
		}

		// Return TRUE if the number of parsed
		// entries is > 0
		return(this->position > 0);
	}

	bool update(bcf_entry_type& entry){
		// Check again because we might use it
		// iteratively at some point
		// i.e. not operating through the
		// build() function
		if(!entry.isBiallelic())
			return false;

		// Cycle over genotypes at this position
		// Ignore phasing at this stage
		//
		// Genotypes are thus:
		// 0/0 -> 0000b = 0
		// 0/1 -> 0001b = 1
		// 0/. -> 0010b = 2
		// 1/0 -> 0100b = 4
		// 1/1 -> 0101b = 5
		// 1/. -> 0110b = 6
		// ./0 -> 1000b = 8
		// ./1 -> 1001b = 9
		// ./. -> 1010b = 10
		//
		// Update GT_array
		U32 j = 0;
		for(U32 i = 0; i < 2*this->n_samples; i += 2, ++j){
			// this->GT_array[j] = PACK(entry.genotypes[i],entry.genotypes[i+1])
		}

		// Build PPA
		// 3^2 = 9 state radix sort over
		// states: alleles \in {00, 01, 11}
		// b entries in a TWK block B
		// This is equivalent to a radix sort
		// on the alphabet {0,1,...,8}
		U32 target_ID = 0;
		for(U32 j = 0; j < this->n_samples; ++j){
			// Determine correct bin
			switch(this->GT_array[this->ppa[j]]){
			case 0:  target_ID = 0; break;
			case 1:  target_ID = 1; break;
			case 2:  target_ID = 2; break;
			case 4:  target_ID = 3; break;
			case 5:  target_ID = 4; break;
			case 6:  target_ID = 5; break;
			case 8:  target_ID = 6; break;
			case 9:  target_ID = 7; break;
			case 10: target_ID = 8; break;
			default: std::cerr << "illegal" << std::endl; exit(1);
			}

			this->bins[target_ID][this->p_i[target_ID]] = this->ppa[j];
			++this->p_i[target_ID];
		} // end loop over individuals at position i

		// Update PPA data
		U32 cum_pos = 0;
		for(U32 i = 0; i < 9; ++i){
			memcpy(&this->ppa[cum_pos], &this->bins[i], this->p_i[i]*sizeof(U32));
			cum_pos += this->p_i[i];
			this->p_i[i] = 0;
		}

		// Keep track of how many entries we've iterated over
		++this->position;

		return true;
	}

	inline const U64& getSamples(void) const{ return(this->n_samples); }
	inline const U32& size(void) const{ return(this->position); }
	inline const U32& operator[](const U32& p) const{return(this->ppa[p]); }

private:
	U64 n_samples;
	U32 position;
	U32* p_i;
	U32* bins;
	//U32 a_i;
	//U32 b_i;
	//U32* a;
	//U32* b;
	U32* ppa;
	BYTE* GT_array; // packed genotype array
};

} /* namespace Algorithm */
} /* namespace Tomahawk */

#endif /* ALGORITHM_COMPRESSION_RADIXSORTGT_H_ */
