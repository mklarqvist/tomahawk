#ifndef GZCONTROLLER_H_
#define GZCONTROLLER_H_

#include <fstream>

#include "../helpers.h"
#include "../third_party/zlib/zlib.h"
#include "BasicBuffer.h"
#include "TGZFHeader.h"

namespace Tomahawk{
namespace IO{

class GZController{
	typedef GZController self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef TGZFHeader header_type;

public:
	GZController();
	GZController(const char* data, const U32 length);
	GZController(const U32 largest_block_size);
	~GZController();

	void Clear();
	bool Inflate(buffer_type& input, buffer_type& output, const header_type& header) const;
	bool Inflate(buffer_type& input, buffer_type& output) const;
	bool __Inflate(buffer_type& input, buffer_type& output, const header_type& header) const;
	bool Deflate(buffer_type& buffer);
	bool Deflate(buffer_type& meta, buffer_type& rle);
	U32 InflateSize(buffer_type& input) const;

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(entry.buffer_.data, entry.buffer_.pointer);
		return stream;
	}

public:
	buffer_type buffer_;
};

}
}

#endif /* GZCONTROLLER_H_ */
