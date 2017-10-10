#ifndef TOTEMPOLEHEADER_H_
#define TOTEMPOLEHEADER_H_

#include <bitset>

namespace Tomahawk {
namespace Totempole {

struct TotempoleHeaderBase{
	typedef TotempoleHeaderBase self_type;

public:
	TotempoleHeaderBase() : version(0), samples(0){}
	TotempoleHeaderBase(const U64 samples) :
						version(Constants::PROGRAM_VERSION),
						samples(samples)
	{}

	friend std::ostream& operator<<(std::ofstream& stream, const self_type& header){
		stream.write(reinterpret_cast<const char*>(&Constants::PROGRAM_VERSION), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&header.samples), sizeof(U64));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		stream.read(reinterpret_cast<char *>(&header.version), sizeof(float));
		stream.read(reinterpret_cast<char *>(&header.samples), sizeof(U64));
		return(stream);
	}

public:
	float version;	// version used to write header
	U64 samples;	// number of samples
};

struct TotempoleHeader : public TotempoleHeaderBase{
	typedef TotempoleHeader self_type;

public:
	TotempoleHeader() :
		controller(0),
		blocks(0),
		largest_uncompressed(0)
	{}

	// This ctor is used during construction
	// only possible when sample count is known
	TotempoleHeader(const U64 samples) :
		TotempoleHeaderBase(samples),
		controller(0),
		blocks(0),
		largest_uncompressed(0)
	{}
	~TotempoleHeader(){}

	friend std::ostream& operator<<(std::ofstream& stream, const self_type& header){
		// version | sample count | controller byte | blocks | largest uncompressed
		stream.write(reinterpret_cast<const char*>(&header.version), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&header.samples), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&header.controller), sizeof(BYTE)); // Controller byte
		// At end-of-file, reopen file as in | out | binary and seek to this position and overwrite with the correct position
		stream.write(reinterpret_cast<const char*>(&header.blocks), sizeof(U32)); // Number of blocks in Tomahawk
		stream.write(reinterpret_cast<const char*>(&header.largest_uncompressed), sizeof(U32)); // Size of largest uncompressed block
		return(stream);
	}

	friend std::ostream& operator<<(std::ostream& os, const self_type& block){
		os <<
		"version: " << block.version << '\n' <<
		"samples: " << block.samples << '\n' <<
		"controller: " << std::bitset<8>(block.controller) << '\n' <<
		"blocks: " << block.blocks << '\n' <<
		"largest: " << block.largest_uncompressed;

		return(os);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& block){
		stream.read(reinterpret_cast<char *>(&block.version), sizeof(float));
		stream.read(reinterpret_cast<char *>(&block.samples), sizeof(U64));
		stream.read(reinterpret_cast<char *>(&block.controller), sizeof(BYTE));
		stream.read(reinterpret_cast<char *>(&block.blocks), sizeof(U32));
		stream.read(reinterpret_cast<char *>(&block.largest_uncompressed), sizeof(U32));

		return(stream);
	}

public:
	BYTE controller;	// controller block
	U32 blocks;			// number of blocks in Tomahawk
	U32 largest_uncompressed;	// largest block-size in bytes
};

}
}

#endif /* TOTEMPOLEHEADER_H_ */
