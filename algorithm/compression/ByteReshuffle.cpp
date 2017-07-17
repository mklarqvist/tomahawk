#include <algorithm>
#include "ByteReshuffle.h"

namespace Tomahawk {
namespace Algorithm{

void ByteReshuffle::SetBitWidth(const BYTE width){
	this->byte_width = width;
	if(this->bits){
		switch(width){
		case 1: this->func = &Tomahawk::Algorithm::ByteReshuffle::ReshuffleBits16; this->shufflers = new BitShuffler[8];  break;
		case 2: this->func = &Tomahawk::Algorithm::ByteReshuffle::ReshuffleBits16; this->shufflers = new BitShuffler[16]; break;
		case 4: this->func = &Tomahawk::Algorithm::ByteReshuffle::ReshuffleBits32; this->shufflers = new BitShuffler[32]; break;
		case 8: this->func = &Tomahawk::Algorithm::ByteReshuffle::ReshuffleBits64; this->shufflers = new BitShuffler[64]; break;
		}
		return;
	}

	switch(width){
	case 1: this->func = &Tomahawk::Algorithm::ByteReshuffle::Reshuffle8;  break;
	case 2: this->func = &Tomahawk::Algorithm::ByteReshuffle::Reshuffle16; break;
	case 4: this->func = &Tomahawk::Algorithm::ByteReshuffle::Reshuffle32; break;
	case 8: this->func = &Tomahawk::Algorithm::ByteReshuffle::Reshuffle64; break;
	}
}

void ByteReshuffle::SetShufflers(const IO::BasicBuffer& input, IO::BasicBuffer& output){
	//this->RankBits(input);

	const U32 bit_width = this->byte_width*8;
	//U32 rank_index[bit_width];
	//U32 rank[bit_width];
	//for(U32 i = 0; i < bit_width; ++i)
	//	rank_index[i] = i;

	struct BitRankHelper{
		U32 index;
		U64 setCounts;

		bool operator<(const BitRankHelper& other) const{
			return this->setCounts < other.setCounts;
		}
	};
	BitRankHelper helpers[bit_width];

	for(U32 i = 0; i < bit_width; ++i){
		helpers[i].index = i;
		helpers[i].setCounts = 0;
	}

	for(U32 i = 0; i < input.size(); i += this->byte_width){
		U32 idx = 0;
		for(U32 b = 0; b < this->byte_width; ++b){
			for(U32 j = 0; j < 8; ++j){
				helpers[idx].setCounts += (input[i+b] & (1 << j)) > 0 ? 1 : 0;
				//std::cerr << i << '/' << input.size() << '\t' << idx << '/' << bit_width << '\t' << i+b << '\t' << std::bitset<8>(1 << j) << std::endl;
				++idx;
			}
		}
	}
	std::sort(&helpers[0], &helpers[bit_width]);
	//for(U32 i = 0; i < bit_width; ++i)
	//	std::cerr << i << '\t' << helpers[i].index << '\t' << helpers[i].setCounts << std::endl;

	for(U32 i = 0; i < bit_width; ++i){
		//std::cerr << helpers[i].index << ": " << (input.size()/bit_width) * i << " bit " << helpers[i].index % 8 << std::endl;
		this->shufflers[helpers[i].index].SetDestination(&output[(input.size()/bit_width) * i], helpers[i].index % 8);
	}

	// resize output based on input
	output.resize(input);
}

void ByteReshuffle::RankBits(const IO::BasicBuffer& input){
	BYTE pos = 0;
	memset(this->frequency_bits, 0, 64);
	for(U32 i = 0; i < input.pointer; ++i){
		for(BYTE j = 0; j < 8; ++j, ++pos){
			this->frequency_bits[pos] += (input[i] & (1 << j)) > 0 ? 1 : 0;
		}

		if(pos == this->byte_width*8)
		 pos = 0;
	}
}

void ByteReshuffle::Reshuffle(IO::BasicBuffer& input, IO::BasicBuffer& output){	(*this.*func)(input, output); }

void ByteReshuffle::Reshuffle8(IO::BasicBuffer& input, IO::BasicBuffer& output){
	std::cerr << "here in 8 repack bytes" << std::endl;
	return;
}

void ByteReshuffle::ReshuffleBits8(IO::BasicBuffer& input, IO::BasicBuffer& output){
	std::cerr << "here in 8 repack bits" << std::endl;
	return;
}

void ByteReshuffle::Reshuffle16(IO::BasicBuffer& input, IO::BasicBuffer& output){
	std::cerr << "here in 16 repack bytes" << std::endl;
	char* startA = &output[0];
	char* startB = &output[input.size()/2];

	for(U32 i = 0; i < input.size(); i += 2){
		*startA = input[i];
		*startB = input[i+1];
	}
}

void ByteReshuffle::ReshuffleBits16(IO::BasicBuffer& input, IO::BasicBuffer& output){
	this->SetShufflers(input, output);

	output.resize(input);
	memset(output.data, 0, input.size()*sizeof(char));
	output.pointer = input.size();

	for(U32 i = 0; i < input.size(); i += 2){
		for(U32 j = 8; j < 16; ++j)
			this->shufflers[j] += input[i+1];

		for(U32 j = 0; j < 8; ++j)
			this->shufflers[j] += input[i];
	}

}

void ByteReshuffle::ReshuffleBits32(IO::BasicBuffer& input, IO::BasicBuffer& output){
	this->SetShufflers(input, output);

	output.resize(input);
	memset(output.data, 0, input.size()*sizeof(char));
	output.pointer = input.size();

	for(U32 i = 0; i < input.size(); i += 4){
		for(U32 j = 24; j < 32; ++j)
			this->shufflers[j] += input[i+3];

		for(U32 j = 16; j < 24; ++j)
			this->shufflers[j] += input[i+2];

		for(U32 j = 8; j < 16; ++j)
			this->shufflers[j] += input[i+1];

		for(U32 j = 0; j < 7; ++j)
			this->shufflers[j] += input[i];
	}

	/*
	// rle pack
	char* bwt_packed_rle = new char[input.pointer-empty_byte_count];
	U32 bwt_packed_rle_pointer = 0;

	U32 run_length = 1;
	char* previous = &bwt_packed[0];
	for(U32 i = 1; i < input.pointer-empty_byte_count; ++i){
		if(bwt_packed[i] != *previous){
			if(run_length >= 4){
				//std::cerr << run_length << '|' << (int)*previous << std::endl;
				U32 value = run_length << 8;
				value ^= *previous;
				U32* target = reinterpret_cast<U32*>(&bwt_packed_rle[bwt_packed_rle_pointer]);
				*target = value;
				bwt_packed_rle_pointer += 4;

			} else {
				for(U32 j = 0; j < run_length; ++j){
					//std::cerr << (int)*previous << std::endl;
					bwt_packed_rle[bwt_packed_rle_pointer] = *previous;
					++bwt_packed_rle_pointer;
				}
			}

			previous = &bwt_packed[i];
			run_length = 1;
			continue;
		}
		++run_length;
	}

	std::cerr << "INput: " << input.pointer << "\tUnpacked: " << input.pointer-empty_byte_count << "\tpacked: " << bwt_packed_rle_pointer << "\t" << (double)input.pointer/bwt_packed_rle_pointer*100 << "%" << std::endl;
	std::cout.write(bwt_packed_rle, bwt_packed_rle_pointer);

	delete [] bwt_packed_rle;
	*/
	//std::cout.write(bwt_packed, input.pointer-empty_byte_count);

	//delete [] bwt_packed;
}


void ByteReshuffle::Reshuffle32(IO::BasicBuffer& input, IO::BasicBuffer& output){
	std::cerr << "here in 32 repack bytes" << std::endl;
	/*
	memset(this->reshuffled_data, 0, block.buffer.pointer);
	// Copy first line
	memcpy(this->reshuffled_data, block.buffer.data,Constants::TOMAHAWK_LINE_HEADER_LENGTH);
	U32 controller_index = Constants::TOMAHAWK_LINE_HEADER_LENGTH;

	const U32 offset = block.variants*Constants::TOMAHAWK_LINE_HEADER_LENGTH;
	const U32 max_width = block.buffer.pointer - offset;

	char* startA = &this->reshuffled_data[0+offset];
	char* startB = &this->reshuffled_data[max_width/4*1+offset];
	char* startC = &this->reshuffled_data[max_width/4*2+offset];
	char* startD = &this->reshuffled_data[max_width/4*3+offset];

	const TomahawkLine* base = reinterpret_cast<const TomahawkLine*>(&block.buffer.data[0]);
	//std::cout << base->position << "\t" << base->runs << "\t" << (int)base->ref << "\t" << (int)base->alt << std::endl;
	U32 index = Constants::TOMAHAWK_LINE_HEADER_LENGTH;
	while(index < block.buffer.pointer){
		for(U32 i = 0; i < base->runs; ++i, index+=4, ++startA, ++startB, ++startC, ++startD){
			//const TomahawkLineEntry<U16>* temp = reinterpret_cast<const TomahawkLineEntry<U16>*>(&data[index]);
			//std::cout << temp->length << ':' << temp->alleleA << ',' << temp->alleleB << '\t' << std::bitset<8>(data[index]) << "\t" << std::bitset<8>(data[index+1]) << std::endl;
			*startA = block.buffer.data[index];
			*startB = block.buffer.data[index+2];
			*startC = block.buffer.data[index+3];
			*startD = block.buffer.data[index+4];
			//std::cout << index << "->" << index+1 << "\t" << std::bitset<8>(*startA) << "\t" << std::bitset<8>(*startB) << std::endl;
		}
		base = reinterpret_cast<const TomahawkLine*>(&block.buffer.data[index]);
		memcpy(&this->reshuffled_data[controller_index], &block.buffer.data[index], Constants::TOMAHAWK_LINE_HEADER_LENGTH);
		controller_index += Constants::TOMAHAWK_LINE_HEADER_LENGTH;
		index += Constants::TOMAHAWK_LINE_HEADER_LENGTH;
		//std::cout << base->position << "\t" << base->runs << "\t" << (int)base->ref << "\t" << (int)base->alt << std::endl;
	}
	//exit(1);
	 */
}

void ByteReshuffle::Reshuffle64(IO::BasicBuffer& input, IO::BasicBuffer& output){
	std::cerr << "here in 64 repack bytes" << std::endl;
}

void ByteReshuffle::ReshuffleBits64(IO::BasicBuffer& input, IO::BasicBuffer& output){
	std::cerr << "here in 64 repack bits" << std::endl;
	return;
}

}
} /* namespace Tomahawk */
