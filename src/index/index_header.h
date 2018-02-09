#ifndef TOTEMPOLEHEADER_H_
#define TOTEMPOLEHEADER_H_

#include <bitset>

namespace Tomahawk {
namespace Totempole {

struct IndexHeader{
public:
	typedef IndexHeader self_type;

public:
	IndexHeader() :
		n_blocks(0),
		largest_uncompressed(0)
	{}

	// This ctor is used during construction
	// only possible when sample count is known
	IndexHeader(const U64 samples) :
		n_blocks(0),
		largest_uncompressed(0)
	{}

	~IndexHeader() = default;

	inline const U32 size(void) const{ return(this->n_blocks); }
	inline const U32 getLargestUncompressedBlock(void) const{ return(this->largest_uncompressed); }

	friend std::ostream& operator<<(std::ofstream& stream, const self_type& header){
		stream.write(reinterpret_cast<const char*>(&header.n_blocks), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.largest_uncompressed), sizeof(U32));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& block){
		stream.read(reinterpret_cast<char *>(&block.n_blocks), sizeof(U32));
		stream.read(reinterpret_cast<char *>(&block.largest_uncompressed), sizeof(U32));

		return(stream);
	}

public:
	U32 n_blocks;               // number of blocks in Tomahawk
	U32 largest_uncompressed;   // largest block-size in bytes
};

}
}

#endif /* TOTEMPOLEHEADER_H_ */
