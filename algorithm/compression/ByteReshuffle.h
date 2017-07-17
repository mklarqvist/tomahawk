#ifndef BYTERESHUFFLE_H_
#define BYTERESHUFFLE_H_

#include "../../TypeDefinitions.h"
#include "BitShuffler.h"
#include "../../io/BasicBuffer.h"
#include "../../tomahawk/MagicConstants.h"

namespace Tomahawk {
namespace Algorithm{

class ByteReshuffle {
	typedef void (ByteReshuffle::*shuffler)(IO::BasicBuffer& input, IO::BasicBuffer& output); // Type cast pointer to function

public:
	ByteReshuffle(const bool bits = true) :
		bits(bits),
		byte_width(1),
		func(nullptr),
		shufflers(nullptr)
	{
	}

	~ByteReshuffle(void){
		delete [] this->shufflers;
	}

	void Reshuffle(IO::BasicBuffer& input, IO::BasicBuffer& output);
	void SetBitWidth(const BYTE width);

private:
	void RankBits(const IO::BasicBuffer& input);
	void SetShufflers(const IO::BasicBuffer& input, IO::BasicBuffer& output);

	void Reshuffle8(IO::BasicBuffer& input, IO::BasicBuffer& output);
	void Reshuffle16(IO::BasicBuffer& input, IO::BasicBuffer& output);
	void Reshuffle32(IO::BasicBuffer& input, IO::BasicBuffer& output);
	void Reshuffle64(IO::BasicBuffer& input, IO::BasicBuffer& output);

	void ReshuffleBits8(IO::BasicBuffer& input, IO::BasicBuffer& output);
	void ReshuffleBits16(IO::BasicBuffer& input, IO::BasicBuffer& output);
	void ReshuffleBits32(IO::BasicBuffer& input, IO::BasicBuffer& output);
	void ReshuffleBits64(IO::BasicBuffer& input, IO::BasicBuffer& output);

//private:
public:
	const bool bits;
	BYTE byte_width;
	shuffler func;
	BitShuffler* shufflers;
	U32 frequency_bits[64];
};


}
} /* namespace Tomahawk */

#endif /* BYTERESHUFFLE_H_ */
