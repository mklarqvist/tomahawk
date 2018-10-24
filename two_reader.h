#ifndef TWO_READER_H_
#define TWO_READER_H_

#include <fstream>

#include "buffer.h"
#include "core.h"
#include "index.h"
#include "zstd_codec.h"
#include "vcf_utils.h"

namespace tomahawk {

class twk1_two_iterator {
public:
	twk1_two_iterator() : stream(nullptr){}
	~twk1_two_iterator(){ }

	bool NextBlockRaw();
	bool NextBlock();
	inline const twk1_two_block_t& GetBlock(void) const{ return(this->blk); }

public:
	ZSTDCodec zcodec; // support codec
	twk_buffer_t buf; // support buffer
	twk_oblock_two_t oblk;// block wrapper
	twk1_two_block_t blk; // block
	std::istream* stream; // stream pointer
};

/**<
 * Reader of twk files.
 */
class two_reader {
public:
	two_reader() : buf(nullptr), stream(nullptr){}
	~two_reader(){
		delete stream;
	}

	/**<
	 * Open a target two file. File header, index, and footer will be read
	 * and parsed as part of the opening procedure. If these pass without errors
	 * then return TRUE. Returns FALSE otherwise.
	 * @param file Input target two file.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool Open(std::string file){

		fstream.open(file, std::ios::in|std::ios::binary|std::ios::ate);
		if(!fstream.good()){
			std::cerr << "failed to open: " << file << std::endl;
			return false;
		}
		buf = fstream.rdbuf();
		stream = new std::istream(buf);

		//stream->open(file, std::ios::in|std::ios::binary|std::ios::ate);
		uint64_t filesize = stream->tellg();
		stream->seekg(0);

		// read magic
		char magic[TOMAHAWK_LD_MAGIC_HEADER_LENGTH];
		stream->read(magic, TOMAHAWK_LD_MAGIC_HEADER_LENGTH);
		if(strncmp(magic, TOMAHAWK_LD_MAGIC_HEADER.data(), TOMAHAWK_LD_MAGIC_HEADER_LENGTH) != 0){
			std::cerr << "failed to read two magic" << std::endl;
			return false;
		}

		// Read, decompress, and parse header
		uint64_t buf_size = 0, obuf_size = 0;
		stream->read(reinterpret_cast<char*>(&buf_size), sizeof(uint64_t));
		stream->read(reinterpret_cast<char*>(&obuf_size),sizeof(uint64_t));
		twk_buffer_t obuf(obuf_size);
		twk_buffer_t buf(buf_size);
		stream->read(obuf.data(),obuf_size);
		obuf.n_chars_ = obuf_size;
		//std::cerr << "header=" << buf_size << "," << obuf_size << "/" << buf.capacity() << "/" << obuf.capacity() << std::endl;

		ZSTDCodec zcodec;
		if(zcodec.Decompress(obuf, buf) == false){
			std::cerr << "failed to decompress header" << std::endl;
			return false;
		}
		//std::cerr << "bufs=" << buf.size() << "==" << buf_size << std::endl;
		assert(buf.size() == buf_size);
		buf >> hdr;
		buf.reset(); obuf.reset();
		//std::cerr << "done hdr" << std::endl;

		// Remember seek point to start of data.
		uint64_t data_start = stream->tellg();
		std::cerr << "start of data=" << (data_start) << std::endl;

		// seek to end-of-file
		// seek back to end of file marker and position where index offset is stored
		stream->seekg(filesize - TOMAHAWK_FILE_EOF_LENGTH - sizeof(uint64_t));
		uint64_t offset_start_index = 0;
		stream->read(reinterpret_cast<char*>(&offset_start_index), sizeof(uint64_t));
		//std::cerr << "seek offset=" << offset_start_index << "/" << filesize << std::endl;

		// Seek to start of offst
		stream->seekg(offset_start_index);
		if(stream->good() == false){
			std::cerr << "failed seek" << std::endl;
			return false;
		}
		//std::cerr << "seek good=" << stream->tellg() << "/" << filesize << std::endl;

		// Load index
		uint8_t marker = 0;
		stream->read(reinterpret_cast<char*>(&marker),   sizeof(uint8_t));
		stream->read(reinterpret_cast<char*>(&buf_size), sizeof(uint64_t));
		stream->read(reinterpret_cast<char*>(&obuf_size),sizeof(uint64_t));
		obuf.resize(obuf_size), buf.resize(buf_size);
		//std::cerr << "before read=" << obuf_size << std::endl;
		stream->read(obuf.data(),obuf_size);
		obuf.n_chars_ = obuf_size;
		//std::cerr << "header=" << buf_size << "," << obuf_size << "/" << buf.capacity() << "/" << obuf.capacity() << std::endl;

		if(zcodec.Decompress(obuf, buf) == false){
			std::cerr << "failed to decompress" << std::endl;
			return false;
		}
		buf >> index;

		// Seek back to the beginning of data.
		std::cerr << "current pos=" << (uint64_t)stream->tellg() << std::endl;
		stream->seekg(data_start);
		std::cerr << "current pos after=" << (uint64_t)stream->tellg() << std::endl;

		return(stream->good());
	}

public:
	std::streambuf* buf;
	std::istream* stream;
	std::ifstream fstream;
	io::VcfHeader hdr;
	IndexOutput index;
};

}

#endif /* TWO_READER_H_ */
