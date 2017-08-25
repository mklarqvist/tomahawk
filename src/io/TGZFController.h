#ifndef GZCONTROLLER_H_
#define GZCONTROLLER_H_

#include <fstream>

#include "../support/helpers.h"
#include "../third_party/zlib/zlib.h"
#include "BasicBuffer.h"
#include "GZFHeader.h"

namespace Tomahawk{
namespace IO{

class TGZFController{
	typedef TGZFController self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef TGZFHeader header_type;

public:
	TGZFController();
	TGZFController(const char* data, const U32 length);
	TGZFController(const U32 largest_block_size);
	~TGZFController();

	void Clear();
	bool Inflate(buffer_type& input, buffer_type& output, const header_type& header) const;
	bool Inflate(buffer_type& input, buffer_type& output) const;
	bool Deflate(const buffer_type& buffer);
	bool Deflate(buffer_type& meta, buffer_type& rle);
	U32 InflateSize(buffer_type& input) const;
	bool InflateBlock(std::ifstream& stream, buffer_type& input);

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(entry.buffer.data, entry.buffer.pointer);
		return stream;
	}

private:
	bool __Inflate(buffer_type& input, buffer_type& output, const header_type& header) const;

public:
	buffer_type buffer;
};

}
}

#endif /* GZCONTROLLER_H_ */
