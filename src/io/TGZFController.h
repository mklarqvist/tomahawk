#ifndef GZCONTROLLER_H_
#define GZCONTROLLER_H_

#include <fstream>

#include "../support/helpers.h"
#include "../third_party/zlib/zconf.h"
#include "../third_party/zlib/zlib.h"
#include "BasicBuffer.h"
#include "GZFHeader.h"

namespace Tomahawk{
namespace IO{

class TGZFController{
private:
	typedef TGZFController self_type;

protected:
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
	bool InflateBlock(std::ifstream& stream, buffer_type& input);

	bool Deflate(const buffer_type& buffer);
	bool Deflate(buffer_type& meta, buffer_type& rle);

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(entry.buffer.data, entry.buffer.pointer);
		return stream;
	}

private:
	bool __Inflate(buffer_type& input, buffer_type& output, const header_type& header) const;

public:
	buffer_type buffer;
};

class TGZFControllerStream : public TGZFController{
public:
	TGZFControllerStream();
	~TGZFControllerStream();

	bool InflateOpen(std::ifstream& stream);
	bool Inflate(std::ifstream& stream, const BYTE* output, const U32& avail_out, U32& return_size);
	bool __Inflate(std::ifstream& stream, const BYTE* output, const U32 avail_out, U32& return_size);

private:
	U32 chunk_size;
	U32 avail_in;
	U32 avail_offset;
	U32 bytes_read;
	U32 BSIZE;
	z_stream d_stream;
};

template <class T>
class TGZFEntryIterator : public TGZFControllerStream{
	typedef TGZFControllerStream parent_type;

public:
	TGZFEntryIterator();
	~TGZFEntryIterator();

	bool nextEntry(T*& entry);
	bool Inflate(std::ifstream& stream, const BYTE* output, const U32& avail_out, U32& return_size);

private:
	void reset(void){ this->pointer = 0; this->n_entries = 0; this->entries = nullptr; }

private:
	U32 pointer;
	U32 n_entries;
	// Offset must equal a TGZF boundary
	// no checks are made
	U64 IO_start_offset; // start TGZF block offset
	U64 IO_end_offset; // end TGZF block offset
	const T* entries;
};

//template <class T>
//bool TGZFEntryIterator<T>::Inflate(std::ifstream& stream, const BYTE* output, const U32& avail_out, U32& return_size);

}
}

#endif /* GZCONTROLLER_H_ */
