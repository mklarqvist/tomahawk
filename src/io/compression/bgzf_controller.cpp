#include "bgzf_controller.h"

#include <limits>
#include <fstream>

#include "gz_constants.h"
#include "third_party/zlib/zconf.h"
#include "third_party/zlib/zlib.h"

namespace tomahawk {
namespace io {


BGZFController::BGZFController(){}

BGZFController::BGZFController(const char* data, const U32 length){}

BGZFController::~BGZFController(){ }

void BGZFController::Clear(){ this->buffer.reset(); }

U32 BGZFController::InflateSize(buffer_type& input) const{
	const header_type& header = *reinterpret_cast<const header_type* const>(input.data());
	if(!header.Validate()){
		 std::cerr << helpers::timestamp("ERROR","BGZF") << "Invalid BGZF header" << std::endl;
		 std::cerr << helpers::timestamp("DEBUG","BGZF") << "Output length: " << header.BSIZE << std::endl;
		 std::cerr << helpers::timestamp("DEBUG","BGZF") << std::endl;
		 std::cerr << header << std::endl;
		 exit(1);
	}

	return header.BSIZE;
}

bool BGZFController::Inflate(buffer_type& input, buffer_type& output) const{
	const header_type& header = *reinterpret_cast<const header_type* const>(&input[0]);
	if(!header.Validate()){
		 std::cerr << helpers::timestamp("ERROR","BGZF") << "Invalid BGZF header" << std::endl;
		 std::cerr << helpers::timestamp("DEBUG","BGZF") << "Output length: " << header.BSIZE << std::endl;
		 std::cerr << helpers::timestamp("DEBUG","BGZF") << std::endl;
		 std::cerr << header << std::endl;
		 exit(1);
	}

	return(this->__Inflate(input, output, header));
}

bool BGZFController::Inflate(buffer_type& input, buffer_type& output, const header_type& header) const{
	return(this->__Inflate(input, output, header));
}

bool BGZFController::__Inflate(buffer_type& input, buffer_type& output, const header_type& header) const{
	const U32& uncompressedLength = *reinterpret_cast<const U32*>(&input[input.size() - sizeof(U32)]);
	if(output.size() + uncompressedLength >= output.capacity())
		output.resize((output.size() + uncompressedLength) + 65536);

	//U32* crc = reinterpret_cast<U32*>(&input.data[input.size() - 2*sizeof(U32)]);

	// Bug fix for ZLIB when overflowing an U32
	U64 avail_out = output.capacity() - output.size();
	if(avail_out > std::numeric_limits<U32>::max())
		avail_out = std::numeric_limits<U32>::max();

	z_stream zs;
	zs.zalloc    = NULL;
	zs.zfree     = NULL;
	zs.next_in   = (Bytef*)&input[constants::BGZF_BLOCK_HEADER_LENGTH];
	zs.avail_in  = (header.BSIZE + 1) - 16;
	zs.next_out  = (Bytef*)&output[output.size()];
	zs.avail_out = (U32)avail_out;

	int status = inflateInit2(&zs, constants::GZIP_WINDOW_BITS);

	if(status != Z_OK){
		std::cerr << helpers::timestamp("ERROR","BGZF") << "Zlib inflateInit failed: " << (int)status << std::endl;
		exit(1);
	}

	// decompress
	status = inflate(&zs, Z_FINISH);
	if(status != Z_STREAM_END){
		inflateEnd(&zs);
		std::cerr << helpers::timestamp("ERROR","BGZF") << "Zlib inflateEnd failed: " << (int)status << std::endl;
		exit(1);
	}

	// finalize
	status = inflateEnd(&zs);
	if(status != Z_OK){
		inflateEnd(&zs);
		std::cerr << helpers::timestamp("ERROR","BGZF") << "Zlib inflateFinalize failed: " << (int)status << std::endl;
		exit(1);
	}

	//if(zs.total_out == 0)
	//	std::cerr << helpers::timestamp("LOG", "BGZF") << "Detected empty BGZF block" << std::endl;

	output.n_chars += zs.total_out;

	return(true);
}

bool BGZFController::InflateBlock(std::ifstream& stream, buffer_type& input){
	input.resize(sizeof(header_type));
	stream.read(input.data(), io::constants::BGZF_BLOCK_HEADER_LENGTH);
	const header_type* h = reinterpret_cast<const header_type*>(input.data());
	input.n_chars = io::constants::BGZF_BLOCK_HEADER_LENGTH;
	if(!h->Validate()){
		std::cerr << helpers::timestamp("ERROR", "BCF") << "Failed to validate!" << std::endl;
		std::cerr << *h << std::endl;
		return false;
	}

	input.resize(h->BSIZE + 1); // make sure all data will fit

	// Recast because if buffer is resized then the pointer address is incorrect
	// resulting in segfault
	h = reinterpret_cast<const header_type*>(input.data());

	stream.read(&input[io::constants::BGZF_BLOCK_HEADER_LENGTH], (h->BSIZE + 1) - io::constants::BGZF_BLOCK_HEADER_LENGTH);
	if(!stream.good()){
		std::cerr << helpers::timestamp("ERROR", "BCF") << "Truncated file..." << std::endl;
		return false;
	}

	input.n_chars = h->BSIZE + 1;
	const U32 uncompressed_size = *reinterpret_cast<const U32*>(&input[input.size() - sizeof(U32)]);
	this->buffer.resize(uncompressed_size + 1);
	this->buffer.reset();

	if(!this->Inflate(input, this->buffer)){
		std::cerr << helpers::timestamp("ERROR", "BCF") << "Failed inflate!" << std::endl;
		return false;
	}

	// BGZF EOF marker
	if(this->buffer.size() == 0)
		return false;

	return true;
}

} /* namespace IO */
} /* namespace Tomahawk */
