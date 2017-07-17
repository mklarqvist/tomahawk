#ifndef ALGORITHM_GENOTYPEBITPACKER_H_
#define ALGORITHM_GENOTYPEBITPACKER_H_

namespace Tomahawk{
namespace Algorithm{

#define PACKER_BITS_BYTE	(sizeof(BYTE)*8)

class GenotypeBitPacker{
public:
	GenotypeBitPacker(BYTE* destination, const BYTE stepsize) :
		stepSize(stepsize),
		steps(PACKER_BITS_BYTE/stepsize),
		mask((1 << stepsize) - 1),
		currentOffset(this->steps),
		destination_pointer(0),
		destination(destination)
	{}
	~GenotypeBitPacker(){}

	template <class T>
	void add(const BYTE& value, const T& repeat_times){
		//std::cerr << "Input: " << (int)value << '/' << repeat_times << '\t' << (int)this->currentOffset << '/' << (int)this->steps << '\t' << this->destination_pointer << std::endl;
		// remove until we hit remainder

		if(this->currentOffset == 0){
			this->currentOffset = this->steps;
			++this->destination_pointer;
		}

		//if(this->currentOffset > 4 || this->currentOffset < 0)
		//	exit(1);

		T remainder = repeat_times;
		if(this->currentOffset != this->steps){
			//std::cerr << "offset != steps: " << (int)this->currentOffset << std::endl;
			// If the number of repeats exceeds what is available
			// keep adding until hitting the end
			if(repeat_times >= this->currentOffset){
				for(S32 i = this->currentOffset; i > 0; --i){
					//std::cerr << ((i - 1) * this->stepSize) << std::endl;
					//std::cerr << 1 << '@' << this->destination_pointer  << '\t' << (int)value << '\t' << (i - 1) * this->stepSize << std::endl;
					this->destination[this->destination_pointer] ^= (value & this->mask) << ((i - 1) * this->stepSize);
				}
				remainder -= this->currentOffset;
				this->currentOffset = this->steps;
				++this->destination_pointer;
			} else {
				for(S32 i = this->currentOffset; i > (this->currentOffset - repeat_times); --i){
					//std::cerr << ((i - 1) * this->stepSize) << std::endl;
					//std::cerr << "1a\t" << i << '/' << (this->currentOffset - repeat_times) << std::endl;
					//std::cerr << "1a" << '@' << this->destination_pointer  << '\t' << (int)value << '\t' << (i - 1) * this->stepSize << '\t' << "offset: " << (int)this->currentOffset << std::endl;
					this->destination[this->destination_pointer] ^= (value & this->mask) << ((i - 1) * this->stepSize);
				}

				remainder -= repeat_times;
				this->currentOffset -= repeat_times;
				//std::cerr << "1a now " << std::bitset<8>(this->destination[this->destination_pointer]) << " " << this->currentOffset << std::endl;
			}

			if(remainder == 0)
				return;
		}

		// Number of complete additions
		const T remainder_times_balanced = remainder/this->steps;
		if(remainder_times_balanced > 0){
			BYTE copy_me = 0;
			//std::cerr << "Balanced: " << remainder_times_balanced << std::endl;
			for(S32 i = this->steps; i > 0; --i)
				copy_me ^= (value & this->mask) << ((i-1)*this->stepSize);

			//std::cerr << 2 << '@' << this->destination_pointer << '\t' << (int)value << '\t' << std::bitset<8>(copy_me) << std::endl;
			// Memcpy fails alignment and produces nonsense! Use a loop
			for(U32 i = 0; i < remainder_times_balanced; ++i){
				this->destination[this->destination_pointer] = copy_me;
				++this->destination_pointer;
			}
			//std::cerr << remainder_times_balanced << '\t' << remainder%this->steps << std::endl;
			remainder -= remainder_times_balanced*this->steps;
		}

		//std::cerr << "remainder: " << repeat_times << " -> "  << remainder << std::endl;
		for(U32 i = 0; i < remainder; ++i){
			//std::cerr << 3<< '@' << this->destination_pointer << '\t' << ((this->currentOffset - 1) * this->stepSize) << std::endl;
			this->destination[this->destination_pointer] ^= (value & this->mask) << ((this->currentOffset - 1) * this->stepSize);
			--this->currentOffset;
		}

		//std::cerr << this->destination_pointer << '\t' << this->currentOffset << std::endl;
	}

	U32 size(void) const{ return(this->destination_pointer*this->steps + (this->steps - this->currentOffset)); }
	void setDestination(BYTE* dest){ this->destination = dest; }

	void operator+=(const BYTE* value);
	void operator<<(const BYTE* value);
	void operator++(void);
	void operator--(void);

	const BYTE& operator[](const U32& p) const{ return(this->destination[p]); }

private:
	const BYTE stepSize;		// step size (in bits)
	const BYTE steps;			// number of steps in byte
	const BYTE mask;			// insert mask
	int currentOffset;
	U32 destination_pointer;
	BYTE* destination;
};

}
}

#endif /* ALGORITHM_GENOTYPEBITPACKER_H_ */
