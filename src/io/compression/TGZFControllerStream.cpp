#include <cassert>

#include "../../support/MagicConstants.h" // for SILENT
#include "TGZFControllerStream.h"

namespace Tomahawk{
namespace IO{

TGZFControllerStream::TGZFControllerStream() : STATE(TGZF_STATE::TGZF_INIT), chunk_size(65536), total_out(0), bytes_read(0), BSIZE(0){}
TGZFControllerStream::~TGZFControllerStream(){}

bool TGZFControllerStream::InflateOpen(std::ifstream& stream){
	this->buffer.reset();
	this->buffer.resize(this->chunk_size);
	this->bytes_read = 0;
	stream.read(&this->buffer.data[0], IO::Constants::TGZF_BLOCK_HEADER_LENGTH);
	const header_type* h = reinterpret_cast<const header_type*>(&this->buffer.data[0]);

	if(!h->Validate()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TGZF") << "Failed to validate!" << std::endl;
		std::cerr << *h << std::endl;
		exit(1);
	}

	this->BSIZE = h->BSIZE - Constants::TGZF_BLOCK_HEADER_LENGTH - Constants::TGZF_BLOCK_FOOTER_LENGTH; // data to read
	this->total_out = 0;

	this->d_stream = z_stream();
	this->d_stream.zalloc   = Z_NULL;
	this->d_stream.zfree    = Z_NULL;
	this->d_stream.opaque   = Z_NULL;
	this->d_stream.avail_in = 0;
	this->d_stream.next_in  = Z_NULL;
	this->STATE = TGZF_STATE::TGZF_HEADER;

	int ret = inflateInit2(&this->d_stream, Constants::GZIP_WINDOW_BITS);
	if (ret != Z_OK){
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "Failed inflatinit" << std::endl;
		this->STATE = TGZF_STATE::TGZF_ERROR;
		return ret;
	}

	return true;
}

bool TGZFControllerStream::Inflate(std::ifstream& stream, const BYTE* output, const U32& avail_out, U32& return_size){
	if(this->STATE == TGZF_INIT)
		this->InflateOpen(stream);

	U32 avail_out_inner = 0;
	U32 ret_inner = 0;

	while(this->__Inflate(stream, &output[avail_out_inner], avail_out - avail_out_inner, ret_inner)){
		return_size += ret_inner;
		avail_out_inner += ret_inner;
	}

	if(this->STATE == TGZF_STATE::TGZF_END)
		return_size += ret_inner;

	if(return_size == 0)
		return false;

	return true;
}

bool TGZFControllerStream::__Inflate(std::ifstream& stream, const BYTE* output, const U32 avail_out, U32& return_size){
	// No space in output
	if(avail_out == 0){
		//std::cerr << "RETURN no space" << std::endl;
		return false;
	}

	// No data available in buffer
	// load some more
	if(this->d_stream.avail_in == 0){ // and bytes read < BSIZE
		this->buffer.reset();

		U32 read_amount = this->chunk_size;
		if(this->bytes_read + this->chunk_size > this->BSIZE)
			read_amount = this->BSIZE - this->bytes_read;

		stream.read(&this->buffer.data[0], read_amount);
		size_t total = stream.gcount();
		this->bytes_read += total;

		//std::cerr << "READ: " << total << "\t" << stream.tellg() << std::endl;
		this->d_stream.avail_in = total;
		this->d_stream.next_in  = (Bytef*)&this->buffer.data[0];

		if(total == 0){
			std::cerr << Helpers::timestamp("WARNING","TGZF") << "Nothing read!" << std::endl;
			return false;
		}
		this->buffer.pointer = total;
	}

	const U32 tot_out = this->d_stream.total_out;
	this->d_stream.next_out  = (Bytef*)output;
	this->d_stream.avail_out = avail_out;

	int status = inflate(&this->d_stream, Z_NO_FLUSH);

	assert(status != Z_STREAM_ERROR);

	if(status != Z_OK && status != Z_STREAM_END){
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "inflate failed: " << (int)status << std::endl;
		exit(1);
	}

	return_size = this->d_stream.total_out - tot_out; // bytes inflated
	this->total_out += return_size;

	if(status != Z_STREAM_END){
		this->STATE = TGZF_STATE::TGZF_OK;
		return true;
	}

	// otherwise its final
	status = inflateEnd(&this->d_stream);
	if(status != Z_OK){
		inflateEnd(&this->d_stream);
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "Zlib inflateFinalize failed: " << (int)status << std::endl;
		exit(1);
	}

	if(this->d_stream.total_out == 0){
		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "TGZF") << "Detected empty TGZF block" << std::endl;
	}

	this->STATE = TGZF_STATE::TGZF_END;
	return false;
}

}
}
