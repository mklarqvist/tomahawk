#include "tgzf_controller.h"

#include <io/compression/gz_constants.h>
#include <fstream>
#include <limits>
#include <cassert>

#include "third_party/zlib/zconf.h"
#include "third_party/zlib/zlib.h"

namespace tomahawk {
namespace io {

TGZFController::TGZFController(){}

TGZFController::TGZFController(const char* data, const U32 length){}

TGZFController::TGZFController(const U32 largest_block_size) : buffer(largest_block_size){}

TGZFController::~TGZFController(){ this->buffer.deleteAll(); }

void TGZFController::Clear(){ this->buffer.reset(); }

bool TGZFController::Inflate(buffer_type& input, buffer_type& output) const{
	const header_type& header = *reinterpret_cast<const header_type* const>(&input[0]);
	if(!header.Validate()){
		 std::cerr << helpers::timestamp("ERROR","TGZF") << "Invalid TGZF header" << std::endl;
		 std::cerr << helpers::timestamp("DEBUG","TGZF") << "Output length: " << header.BSIZE << std::endl;
		 std::cerr << helpers::timestamp("DEBUG","TGZF") << std::endl;
		 std::cerr << header << std::endl;
		 exit(1);
	}

	return(this->__Inflate(input, output, header));
}

bool TGZFController::Inflate(buffer_type& input, buffer_type& output, const header_type& header) const{
	return(this->__Inflate(input, output, header));
}

bool TGZFController::__Inflate(buffer_type& input, buffer_type& output, const header_type& header) const{
	const U32& uncompressedLength = *reinterpret_cast<const U32*>(&input[input.size() - sizeof(U32)]);
	if(output.size() + uncompressedLength >= output.capacity())
		output.resize((output.size() + uncompressedLength) + 65536);

	// Not used
	//U32* crc = reinterpret_cast<U32*>(&input.data[input.size() - 2*sizeof(U32)]);

	// Fix for ZLIB when overflowing an U32
	U64 avail_out = output.capacity() - output.size();
	if(avail_out > std::numeric_limits<U32>::max())
		avail_out = std::numeric_limits<U32>::max();

	z_stream zs;
	zs.zalloc    = NULL;
	zs.zfree     = NULL;
	zs.next_in   = (Bytef*)&input[constants::TGZF_BLOCK_HEADER_LENGTH];
	zs.avail_in  = (header.BSIZE + 1) - 16;
	zs.next_out  = (Bytef*)&output[output.size()];
	zs.avail_out = (U32)avail_out;

	int status = inflateInit2(&zs, constants::GZIP_WINDOW_BITS);

	if(status != Z_OK){
		std::cerr << helpers::timestamp("ERROR","TGZF") << "Zlib inflateInit failed: " << (int)status << std::endl;
		exit(1);
	}

	// decompress
	status = inflate(&zs, Z_FINISH);
	if(status != Z_STREAM_END){
		inflateEnd(&zs);
		std::cerr << helpers::timestamp("ERROR","TGZF") << "Zlib inflateEnd failed: " << (int)status << std::endl;
		exit(1);
	}

	// finalize
	status = inflateEnd(&zs);
	if(status != Z_OK){
		inflateEnd(&zs);
		std::cerr << helpers::timestamp("ERROR","TGZF") << "Zlib inflateFinalize failed: " << (int)status << std::endl;
		exit(1);
	}

	if(zs.total_out == 0)
		std::cerr << helpers::timestamp("LOG", "TGZF") << "Detected empty TGZF block" << std::endl;

	output.n_chars += zs.total_out;

	return(true);
}

bool TGZFController::Deflate(const buffer_type& buffer){
	this->buffer.resize(buffer);

	memset(this->buffer.data(), 0, constants::TGZF_BLOCK_HEADER_LENGTH);

	this->buffer[0]  = constants::GZIP_ID1;
	this->buffer[1]  = constants::GZIP_ID2;
	this->buffer[2]  = constants::CM_DEFLATE;
	this->buffer[3]  = constants::FLG_FEXTRA;
	this->buffer[9]  = constants::OS_UNKNOWN;
	this->buffer[10] = constants::TGZF_XLEN;
	this->buffer[12] = constants::TGZF_ID1;
	this->buffer[13] = constants::TGZF_ID2;
	this->buffer[14] = constants::TGZF_LEN;
	//buffer 16->20 is set below

	// set compression level
	const int compressionLevel = Z_DEFAULT_COMPRESSION;
	//const int compressionLevel = 9;

	// initialize zstream values
    z_stream zs;
    zs.zalloc    = NULL;
    zs.zfree     = NULL;
    zs.next_in   = (Bytef*)buffer.data();
    zs.avail_in  = buffer.size();
    zs.next_out  = (Bytef*)&this->buffer[constants::TGZF_BLOCK_HEADER_LENGTH];
    zs.avail_out = this->buffer.width -
                   constants::TGZF_BLOCK_HEADER_LENGTH -
                   constants::TGZF_BLOCK_FOOTER_LENGTH;

	// Initialise the zlib compression algorithm
	int status = deflateInit2(&zs,
							  compressionLevel,
							  Z_DEFLATED,
							  constants::GZIP_WINDOW_BITS,
							  constants::Z_DEFAULT_MEM_LEVEL,
							  Z_DEFAULT_STRATEGY);

	if ( status != Z_OK ){
		std::cerr << helpers::timestamp("ERROR", "ZLIB") << "DeflateBlock: zlib deflateInit2 failed" << std::endl;
		return false;
	}

	// compress the data
	status = deflate(&zs, Z_FINISH);

	// if not at stream end
	if ( status != Z_STREAM_END ) {
		deflateEnd(&zs);

		// there was not enough space available in buffer
		std::cerr << helpers::timestamp("ERROR", "ZLIB") << "DeflateBlock: zlib deflate failed (insufficient space)" << std::endl;
		return false;
	}

	// finalize the compression routine
	status = deflateEnd(&zs);
	if ( status != Z_OK ){
		std::cerr << helpers::timestamp("ERROR", "ZLIB") << "DeflateBlock: zlib deflateEnd failed (not ok)" << std::endl;
		return false;
	}

	// update compressedLength
	const U32 compressedLength = zs.total_out +
					       	     constants::TGZF_BLOCK_HEADER_LENGTH +
								 constants::TGZF_BLOCK_FOOTER_LENGTH;

	// store the compressed length
	U32* test = reinterpret_cast<U32*>(&this->buffer[16]);
	*test = compressedLength;
	//std::cerr << helpers::timestamp("DEBUG") << data.pointer << "->" << compressedLength-1 << " stored: " << *test << std::endl;

	std::time_t result = std::time(nullptr);
	std::asctime(std::localtime(&result));
	U32* time = reinterpret_cast<U32*>(&this->buffer[4]);
	*time = result;
	//std::cerr << helpers::timestamp("DEBUG") << "Time: " << *time << std::endl;


	memset(&buffer.buffer[compressedLength - constants::TGZF_BLOCK_FOOTER_LENGTH], 0, constants::TGZF_BLOCK_FOOTER_LENGTH);

	// store the CRC32 checksum
	U32 crc = crc32(0, NULL, 0);
	crc = crc32(crc, (Bytef*)buffer.data(), buffer.size());
	U32* c = reinterpret_cast<U32*>(&this->buffer[compressedLength - constants::TGZF_BLOCK_FOOTER_LENGTH]);
	*c = crc;
	U32 convert = buffer.size(); // avoid potential problems when casting from U64 to U32 by interpretation
	U32* uncompressed = reinterpret_cast<U32*>(&this->buffer[compressedLength - sizeof(U32)]);
	*uncompressed = convert; // Store uncompressed length

	this->buffer.n_chars = compressedLength;
	//std::cerr << "Writing: " << convert << '/' << *uncompressed << '\t' << compressedLength << '\t' << *test << '\t' << buffer.size() << '\t' << "At pos: " << (compressedLength - sizeof(U32)) << '\t' << buffer.pointer << '\t' << *c << '\t' << convert << std::endl;

	return true;
}

bool TGZFController::Deflate(buffer_type& meta, buffer_type& rle){
	meta += rle;
	return(this->Deflate(meta));
}

bool TGZFController::InflateBlock(std::istream& stream, buffer_type& input){
	input.resize(sizeof(header_type));
	stream.read(input.data(), io::constants::TGZF_BLOCK_HEADER_LENGTH);
	const header_type* h = reinterpret_cast<const header_type*>(input.data());
	input.n_chars = io::constants::TGZF_BLOCK_HEADER_LENGTH;
	if(!h->Validate()){
		std::cerr << helpers::timestamp("ERROR", "TGZF") << "Failed to validate!" << std::endl;
		std::cerr << *h << std::endl;
		return false;
	}

	input.resize(h->BSIZE); // make sure all data will fit

	// Recast because if buffer is resized then the pointer address is incorrect
	// resulting in segfault
	h = reinterpret_cast<const header_type*>(input.data());

	stream.read(&input[io::constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - io::constants::TGZF_BLOCK_HEADER_LENGTH);
	if(!stream.good()){
		std::cerr << helpers::timestamp("ERROR", "TGZF") << "Truncated file..." << std::endl;
		return false;
	}

	input.n_chars = h->BSIZE;
	const U32 uncompressed_size = *reinterpret_cast<const U32*>(&input[input.size() - sizeof(U32)]);
	this->buffer.resize(uncompressed_size);
	this->buffer.reset();

	if(!this->Inflate(input, this->buffer)){
		std::cerr << helpers::timestamp("ERROR", "TGZF") << "Failed inflate!" << std::endl;
		return false;
	}

	// TGZF EOF marker
	if(this->buffer.size() == 0)
		return false;

	return true;
}


}
}
