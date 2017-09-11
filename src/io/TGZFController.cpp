#include <fstream>
#include <limits>
#include <cassert>

#include "../third_party/zlib/zconf.h"
#include "../third_party/zlib/zlib.h"
#include "IOConstants.h"
#include "TGZFController.h"

namespace Tomahawk {
namespace IO {

TGZFController::TGZFController(){}

TGZFController::TGZFController(const char* data, const U32 length){}

TGZFController::TGZFController(const U32 largest_block_size) : buffer(largest_block_size){}

TGZFController::~TGZFController(){ this->buffer.deleteAll(); }

void TGZFController::Clear(){ this->buffer.reset(); }

bool TGZFController::Inflate(buffer_type& input, buffer_type& output) const{
	const header_type& header = *reinterpret_cast<const header_type* const>(&input[0]);
	if(!header.Validate()){
		 std::cerr << Helpers::timestamp("ERROR","TGZF") << "Invalid TGZF header" << std::endl;
		 std::cerr << Helpers::timestamp("DEBUG","TGZF") << "Output length: " << header.BSIZE << std::endl;
		 std::cerr << Helpers::timestamp("DEBUG","TGZF") << std::endl;
		 std::cerr << header << std::endl;
		 exit(1);
	}

	return(this->__Inflate(input, output, header));
}

bool TGZFController::Inflate(buffer_type& input, buffer_type& output, const header_type& header) const{
	return(this->__Inflate(input, output, header));
}

bool TGZFController::__Inflate(buffer_type& input, buffer_type& output, const header_type& header) const{
	const U32& uncompressedLength = *reinterpret_cast<const U32*>(&input.data[input.size() - sizeof(U32)]);
	if(output.size() + uncompressedLength >= output.capacity())
		output.resize((output.size() + uncompressedLength) + 65536);

	//U32* crc = reinterpret_cast<U32*>(&input.data[input.size() - 2*sizeof(U32)]);

	// Fix for ZLIB when overflowing an U32
	U64 avail_out = output.capacity() - output.size();
	if(avail_out > std::numeric_limits<U32>::max())
		avail_out = std::numeric_limits<U32>::max();

	z_stream zs;
	zs.zalloc    = NULL;
	zs.zfree     = NULL;
	zs.next_in   = (Bytef*)&input.data[Constants::TGZF_BLOCK_HEADER_LENGTH];
	zs.avail_in  = (header.BSIZE + 1) - 16;
	zs.next_out  = (Bytef*)&output.data[output.pointer];
	zs.avail_out = (U32)avail_out;

	int status = inflateInit2(&zs, Constants::GZIP_WINDOW_BITS);

	if(status != Z_OK){
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "Zlib inflateInit failed: " << (int)status << std::endl;
		exit(1);
	}

	// decompress
	status = inflate(&zs, Z_FINISH);
	if(status != Z_STREAM_END){
		inflateEnd(&zs);
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "Zlib inflateEnd failed: " << (int)status << std::endl;
		exit(1);
	}

	// finalize
	status = inflateEnd(&zs);
	if(status != Z_OK){
		inflateEnd(&zs);
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "Zlib inflateFinalize failed: " << (int)status << std::endl;
		exit(1);
	}

	if(zs.total_out == 0)
		std::cerr << Helpers::timestamp("LOG", "TGZF") << "Detected empty TGZF block" << std::endl;

	output.pointer += zs.total_out;

	return(true);
}

bool TGZFController::Deflate(const buffer_type& buffer){
	this->buffer.resize(buffer);

	memset(this->buffer.data, 0, Constants::TGZF_BLOCK_HEADER_LENGTH);

	this->buffer[0]  = Constants::GZIP_ID1;
	this->buffer[1]  = Constants::GZIP_ID2;
	this->buffer[2]  = Constants::CM_DEFLATE;
	this->buffer[3]  = Constants::FLG_FEXTRA;
	this->buffer[9]  = Constants::OS_UNKNOWN;
	this->buffer[10] = Constants::TGZF_XLEN;
	this->buffer[12] = Constants::TGZF_ID1;
	this->buffer[13] = Constants::TGZF_ID2;
	this->buffer[14] = Constants::TGZF_LEN;
	//buffer 16->20 is set below

	// set compression level
	const int compressionLevel = Z_DEFAULT_COMPRESSION;
	//const int compressionLevel = 9;

	// initialize zstream values
    z_stream zs;
    zs.zalloc    = NULL;
    zs.zfree     = NULL;
    zs.next_in   = (Bytef*)buffer.data;
    zs.avail_in  = buffer.pointer;
    zs.next_out  = (Bytef*)&this->buffer[Constants::TGZF_BLOCK_HEADER_LENGTH];
    zs.avail_out = this->buffer.width -
                   Constants::TGZF_BLOCK_HEADER_LENGTH -
                   Constants::TGZF_BLOCK_FOOTER_LENGTH;

	// Initialise the zlib compression algorithm
	int status = deflateInit2(&zs,
							  compressionLevel,
							  Z_DEFLATED,
							  Constants::GZIP_WINDOW_BITS,
							  Constants::Z_DEFAULT_MEM_LEVEL,
							  Z_DEFAULT_STRATEGY);

	if ( status != Z_OK ){
		std::cerr << Helpers::timestamp("ERROR", "ZLIB") << "DeflateBlock: zlib deflateInit2 failed" << std::endl;
		return false;
	}

	// compress the data
	status = deflate(&zs, Z_FINISH);

	// if not at stream end
	if ( status != Z_STREAM_END ) {
		deflateEnd(&zs);

		// there was not enough space available in buffer
		std::cerr << Helpers::timestamp("ERROR", "ZLIB") << "DeflateBlock: zlib deflate failed (insufficient space)" << std::endl;
		return false;
	}

	// finalize the compression routine
	status = deflateEnd(&zs);
	if ( status != Z_OK ){
		std::cerr << Helpers::timestamp("ERROR", "ZLIB") << "DeflateBlock: zlib deflateEnd failed (not ok)" << std::endl;
		return false;
	}

	// update compressedLength
	const U32 compressedLength = zs.total_out +
					       	     Constants::TGZF_BLOCK_HEADER_LENGTH +
								 Constants::TGZF_BLOCK_FOOTER_LENGTH;

	// store the compressed length
	U32* test = reinterpret_cast<U32*>(&this->buffer[16]);
	*test = compressedLength;
	//std::cerr << Helpers::timestamp("DEBUG") << data.pointer << "->" << compressedLength-1 << " stored: " << *test << std::endl;

	std::time_t result = std::time(nullptr);
	std::asctime(std::localtime(&result));
	U32* time = reinterpret_cast<U32*>(&this->buffer[4]);
	*time = result;
	//std::cerr << Helpers::timestamp("DEBUG") << "Time: " << *time << std::endl;


	memset(&buffer.data[compressedLength - Constants::TGZF_BLOCK_FOOTER_LENGTH], 0, Constants::TGZF_BLOCK_FOOTER_LENGTH);

	// store the CRC32 checksum
	U32 crc = crc32(0, NULL, 0);
	crc = crc32(crc, (Bytef*)buffer.data, buffer.pointer);
	U32* c = reinterpret_cast<U32*>(&this->buffer[compressedLength - Constants::TGZF_BLOCK_FOOTER_LENGTH]);
	*c = crc;
	U32 convert = buffer.pointer; // avoid potential problems when casting from U64 to U32 by interpretation
	U32* uncompressed = reinterpret_cast<U32*>(&this->buffer[compressedLength - sizeof(U32)]);
	*uncompressed = convert; // Store uncompressed length

	this->buffer.pointer = compressedLength;
	//std::cerr << "Writing: " << convert << '/' << *uncompressed << '\t' << compressedLength << '\t' << *test << '\t' << buffer.size() << '\t' << "At pos: " << (compressedLength - sizeof(U32)) << '\t' << buffer.pointer << '\t' << *c << '\t' << convert << std::endl;

	return true;
}

bool TGZFController::Deflate(buffer_type& meta, buffer_type& rle){
	meta += rle;
	return(this->Deflate(meta));
}

bool TGZFController::InflateBlock(std::ifstream& stream, buffer_type& input){
	input.resize(sizeof(header_type));
	stream.read(&input.data[0], IO::Constants::TGZF_BLOCK_HEADER_LENGTH);
	const header_type* h = reinterpret_cast<const header_type*>(&input.data[0]);
	input.pointer = IO::Constants::TGZF_BLOCK_HEADER_LENGTH;
	if(!h->Validate()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TGZF") << "Failed to validate!" << std::endl;
		std::cerr << *h << std::endl;
		return false;
	}

	input.resize(h->BSIZE); // make sure all data will fit

	// Recast because if buffer is resized then the pointer address is incorrect
	// resulting in segfault
	h = reinterpret_cast<const header_type*>(&input.data[0]);

	stream.read(&input.data[IO::Constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - IO::Constants::TGZF_BLOCK_HEADER_LENGTH);
	if(!stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TGZF") << "Truncated file..." << std::endl;
		return false;
	}

	input.pointer = h->BSIZE;
	const U32 uncompressed_size = *reinterpret_cast<const U32*>(&input[input.pointer -  sizeof(U32)]);
	this->buffer.resize(uncompressed_size);
	this->buffer.reset();

	if(!this->Inflate(input, this->buffer)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TGZF") << "Failed inflate!" << std::endl;
		return false;
	}

	// TGZF EOF marker
	if(this->buffer.size() == 0)
		return false;

	return true;
}

///

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
		std::cerr << "failed inflatinit" << std::endl;
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

	if(this->d_stream.total_out == 0)
		std::cerr << Helpers::timestamp("LOG", "TGZF") << "Detected empty TGZF block" << std::endl;

	this->STATE = TGZF_STATE::TGZF_END;
	return false;
}


}
}
