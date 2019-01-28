#include "twk_reader.h"

namespace tomahawk {

/****************************
*  twk1_blk_iterator
****************************/
bool twk1_blk_iterator::NextBlockRaw(){
	if(stream->good() == false){
		std::cerr << "stream died" << std::endl;
		return false;
	}

	uint8_t marker = 0;
	DeserializePrimitive(marker, *stream);
	if(marker == 0){
		//std::cerr << "0 marker found. stopping" << std::endl;
		return false;
	}
	assert(marker == 1);

	*stream >> oblk;
	if(stream->good() == false){
		std::cerr << "stream died" << std::endl;
		return false;
	}

	assert(oblk.bytes.size() == oblk.nc);
	buf.resize(oblk.n);

	return true;
}

bool twk1_blk_iterator::NextBlock(){
	if(this->NextBlockRaw() == false)
		return false;

	// Decompress data
	zcodec.Decompress(oblk.bytes, buf);
	buf >> blk;
	buf.reset();

	return true;
}

/****************************
*  twk_reader
****************************/
bool twk_reader::Open(std::string file){
	fstream.open(file, std::ios::in|std::ios::binary|std::ios::ate);
	if(!fstream.good()){
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to open \"" << file << "\"!" << std::endl;
		return false;
	}
	buf = fstream.rdbuf();
	stream = new std::istream(buf);

	//stream->open(file, std::ios::in|std::ios::binary|std::ios::ate);
	uint64_t filesize = stream->tellg();
	stream->seekg(0);

	// read magic
	char magic[TOMAHAWK_MAGIC_HEADER_LENGTH];
	stream->read(magic, TOMAHAWK_MAGIC_HEADER_LENGTH);
	if(strncmp(magic, TOMAHAWK_MAGIC_HEADER.data(), TOMAHAWK_MAGIC_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to read MAGIC!" << std::endl;
		return false;
	}

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
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to decompress header!" << std::endl;
		return false;
	}
	//std::cerr << "bufs=" << buf.size() << "==" << buf_size << std::endl;
	assert(buf.size() == buf_size);
	buf >> hdr;
	buf.reset(); obuf.reset();
	//std::cerr << "done hdr" << std::endl;

	uint64_t data_start = stream->tellg();

	// seek to end-of-file
	stream->seekg(filesize - TOMAHAWK_FILE_EOF_LENGTH - sizeof(uint64_t));
	uint64_t offset_start_index = 0;
	stream->read(reinterpret_cast<char*>(&offset_start_index), sizeof(uint64_t));
	//std::cerr << "seek offset=" << offset_start_index << "/" << filesize << std::endl;
	stream->seekg(offset_start_index);
	if(stream->good() == false){
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to seek in file!" << std::endl;
		return false;
	}
	//std::cerr << "seek good=" << stream->tellg() << "/" << filesize << std::endl;

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
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to decompress index!" << std::endl;
		return false;
	}
	buf >> index;

	// Seek back to start of data.
	stream->seekg(data_start);

	return true;
}

}
