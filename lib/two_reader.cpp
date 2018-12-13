#include "two_reader.h"

namespace tomahawk {

bool twk1_two_iterator::NextBlockRaw(){
	if(stream->good() == false){
		std::cerr << "stream died" << std::endl;
		return false;
	}

	//std::cerr << "pos=" << stream->tellg() << std::endl;

	uint8_t marker = 0;
	DeserializePrimitive(marker, *stream);
	if(marker == 0){
		//std::cerr << "0 marker found. stopping" << std::endl;
		return false;
	}
	//std::cerr << "marker=" << (int)marker << std::endl;
	//assert(marker == 1);
	if(marker != 1){
		std::cerr << "Marker!=1 is " << (int)marker << " @ " << stream->tellg() << " good=" << stream->good() << std::endl;
		exit(1);
	}

	*stream >> oblk;
	if(stream->good() == false){
		std::cerr << "stream died" << std::endl;
		return false;
	}

	assert(oblk.bytes.size() == oblk.nc);
	buf.resize(oblk.n);
	offset = 0;
	rcd = nullptr;

	return true;
}

bool twk1_two_iterator::NextBlock(){
	if(this->NextBlockRaw() == false)
		return false;

	// Decompress data
	zcodec.Decompress(oblk.bytes, buf);
	buf >> blk;
	buf.reset();
	if(blk.n) rcd = &blk.rcds[0];

	return true;
}

bool twk1_two_iterator::NextRecord(){
	if(offset == blk.n){
		if(this->NextBlock() == false)
			return false;

		offset = 0;
	}
	rcd = &blk.rcds[offset++];
	return true;
}

bool two_reader::Open(std::string file){
	fstream.open(file, std::ios::in|std::ios::binary|std::ios::ate);
	if(!fstream.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to open: " << file << std::endl;
		return false;
	}
	buf = fstream.rdbuf();
	stream = new std::istream(buf);

	uint64_t filesize = stream->tellg();
	stream->seekg(0);

	// read magic
	char magic[TOMAHAWK_LD_MAGIC_HEADER_LENGTH];
	stream->read(magic, TOMAHAWK_LD_MAGIC_HEADER_LENGTH);
	if(strncmp(magic, TOMAHAWK_LD_MAGIC_HEADER.data(), TOMAHAWK_LD_MAGIC_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR") << "Failed to read TWO magic string!" << std::endl;
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

	if(zcodec.Decompress(obuf, buf) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to decompress header!" << std::endl;
		return false;
	}
	assert(buf.size() == buf_size);
	buf >> hdr;
	buf.reset(); obuf.reset();

	// Remember seek point to start of data.
	uint64_t data_start = stream->tellg();

	// seek to end-of-file
	// seek back to end of file marker and position where index offset is stored
	stream->seekg(filesize - TOMAHAWK_FILE_EOF_LENGTH - sizeof(uint64_t));
	uint64_t offset_start_index = 0;
	stream->read(reinterpret_cast<char*>(&offset_start_index), sizeof(uint64_t));

	// Seek to start of offst
	stream->seekg(offset_start_index);
	if(stream->good() == false){
		std::cerr << utility::timestamp("ERROR") << "Failed seek in file!" << std::endl;
		return false;
	}

	// Load index
	uint8_t marker = 0;
	stream->read(reinterpret_cast<char*>(&marker),   sizeof(uint8_t));
	stream->read(reinterpret_cast<char*>(&buf_size), sizeof(uint64_t));
	stream->read(reinterpret_cast<char*>(&obuf_size),sizeof(uint64_t));
	obuf.resize(obuf_size), buf.resize(buf_size);
	stream->read(obuf.data(),obuf_size);
	obuf.n_chars_ = obuf_size;

	if(zcodec.Decompress(obuf, buf) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to decompress!" << std::endl;
		return false;
	}
	buf >> index;

	// Seek back to the beginning of data.
	stream->seekg(data_start);

	// Assign stream.
	it.stream = stream;

	return(stream->good());
}

}
