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
		std::cerr << "0 marker found. stopping" << std::endl;
		return false;
	}
	//std::cerr << "marker=" << (int)marker << std::endl;
	assert(marker == 1);

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

}
