#ifndef BITSHUFFLER_H_
#define BITSHUFFLER_H_

#include "../../TypeDefinitions.h"

namespace Tomahawk {
namespace Algorithm {

class BitShuffler{
public:
	BitShuffler(void) : offsetDistance(1), offsetMask(1), offset(0), pointer(0), destination(nullptr){}
	BitShuffler(char* destination, const BYTE offsetDistance) :
		offsetDistance(offsetDistance),
		offsetMask(1 << this->offsetDistance),
		offset(0),
		pointer(0),
		destination(destination)
	{}
	~BitShuffler(){}

	void operator+=(const char& input){
		if(this->offset == 8){
			this->offset = 0;
			//std::cout << this->pointer << "\tfinal: " << std::bitset<8>(this->destination[this->pointer]) << std::endl;
			++this->pointer;
		}
		this->destination[this->pointer] <<= 1;
		this->destination[this->pointer] ^= (input & this->offsetMask) >> this->offsetDistance;
		//std::cout << this->pointer << "\t" << std::bitset<8>(this->offsetMask) << "\t"  << (int)this->offsetDistance << "\t" << std::bitset<8>((input & this->offsetMask) >> this->offsetDistance) << '\t' << std::bitset<8>(input) << std::endl;

		++this->offset;
	}

	void SetDestination(char* destination, const BYTE offsetDistance){
		this->destination = destination;
		this->offsetDistance = offsetDistance;
		this->offsetMask = 1 << this->offsetDistance;
		this->pointer = 0; // Reset
	}

private:
	BYTE offsetDistance;
	BYTE offsetMask;
	BYTE offset;
	U32 pointer;
	char* destination;
};

}
} /* namespace Tomahawk */

#endif /* BITSHUFFLER_H_ */
