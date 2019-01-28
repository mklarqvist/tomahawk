#ifndef TWK_GENOTYPE_BITPACKER_H_
#define TWK_GENOTYPE_BITPACKER_H_

#include <cstdint>

namespace tomahawk {

/**<
 * This class packs run-length encoded genotypes
 * into a 2-bit encoding system used in vectorized
 * functions throughout Tomahawk
 */
class GenotypeBitPacker {
public:
	const static uint32_t step_size = 2;
	const static uint32_t steps     = 8 / step_size;
	const static uint32_t mask      = ((1<<step_size)-1);

public:
	GenotypeBitPacker(uint8_t* destination, const uint8_t stepsize = 2) :
		currentOffset(steps),
		n(0),
		dst(destination)
	{}

	~GenotypeBitPacker(){}

	void Add(const uint8_t& value, const uint32_t repeat_times){
		if(this->currentOffset == 0){
			this->currentOffset = this->steps;
			++this->n;
		}

		uint32_t remainder = repeat_times;
		if(this->currentOffset != this->steps){
			// If the number of repeats exceeds what is available
			// keep adding until hitting the end
			if(repeat_times >= this->currentOffset){
				for(int i = this->currentOffset; i > 0; --i)
					this->dst[this->n] ^= (value & this->mask) << ((i - 1) * step_size);

				remainder -= this->currentOffset;
				this->currentOffset = this->steps;
				++this->n;
			} else {
				for(int i = this->currentOffset; i > (this->currentOffset - repeat_times); --i)
					this->dst[this->n] ^= (value & this->mask) << ((i - 1) * step_size);


				remainder -= repeat_times;
				this->currentOffset -= repeat_times;
			}

			if(remainder == 0)
				return;
		}

		// Number of complete additions
		const uint32_t remainder_times_balanced = remainder/this->steps;
		if(remainder_times_balanced > 0){
			uint8_t copy_me = 0;
			for(int i = this->steps; i > 0; --i)
				copy_me ^= (value & this->mask) << ((i-1)*step_size);

			// Memcpy fails alignment and produces nonsense! Use a loop
			for(uint32_t i = 0; i < remainder_times_balanced; ++i){
				this->dst[this->n] = copy_me;
				++this->n;
			}
			remainder -= remainder_times_balanced*this->steps;
		}

		for(uint32_t i = 0; i < remainder; ++i){
			this->dst[this->n] ^= (value & this->mask) << ((this->currentOffset - 1) * step_size);
			--this->currentOffset;
		}
	}

	inline uint32_t size(void) const{ return(this->n*this->steps + (this->steps - this->currentOffset)); }
	inline void SetDestination(uint8_t* dest){ this->dst = dest; }
	inline const uint8_t& operator[](const uint32_t& p) const{ return(this->dst[p]); }

public:
	int currentOffset;
	uint32_t n;
	uint8_t* dst;
};

}


#endif /* GENOTYPE_BITPACKER_H_ */
