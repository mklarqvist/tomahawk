#include <limits>

#include "../third_party/zlib/zlib.h"
#include "IOConstants.h"
#include "BGZFController.h"

namespace Tomahawk {
namespace IO {


BGZFController::BGZFController(){}

BGZFController::BGZFController(const char* data, const U32 length){}

BGZFController::~BGZFController(){ this->buffer.deleteAll(); }

void BGZFController::Clear(){ this->buffer.reset(); }

U32 BGZFController::InflateSize(buffer_type& input) const{
	const header_type& header = *reinterpret_cast<const header_type* const>(&input.data[0]);
	if(!header.Validate()){
		 std::cerr << Helpers::timestamp("ERROR","BGZF") << "Invalid BGZF header" << std::endl;
		 std::cerr << Helpers::timestamp("DEBUG","BGZF") << "Output length: " << header.BSIZE << std::endl;
		 std::cerr << Helpers::timestamp("DEBUG","BGZF") << std::endl;
		 std::cerr << header << std::endl;
		 exit(1);
	}

	return header.BSIZE;
}

bool BGZFController::Inflate(buffer_type& input, buffer_type& output) const{
	const header_type& header = *reinterpret_cast<const header_type* const>(&input[0]);
	if(!header.Validate()){
		 std::cerr << Helpers::timestamp("ERROR","BGZF") << "Invalid BGZF header" << std::endl;
		 std::cerr << Helpers::timestamp("DEBUG","BGZF") << "Output length: " << header.BSIZE << std::endl;
		 std::cerr << Helpers::timestamp("DEBUG","BGZF") << std::endl;
		 std::cerr << header << std::endl;
		 exit(1);
	}

	return(this->__Inflate(input, output, header));
}

bool BGZFController::Inflate(buffer_type& input, buffer_type& output, const header_type& header) const{
	return(this->__Inflate(input, output, header));
}

bool BGZFController::__Inflate(buffer_type& input, buffer_type& output, const header_type& header) const{
	const U32& uncompressedLength = *reinterpret_cast<const U32*>(&input.data[input.size() - sizeof(U32)]);
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
	zs.next_in   = (Bytef*)&input.data[Constants::BGZF_BLOCK_HEADER_LENGTH];
	zs.avail_in  = (header.BSIZE + 1) - 16;
	zs.next_out  = (Bytef*)&output.data[output.pointer];
	zs.avail_out = (U32)avail_out;

	int status = inflateInit2(&zs, Constants::GZIP_WINDOW_BITS);

	if(status != Z_OK){
		std::cerr << Helpers::timestamp("ERROR","BGZF") << "Zlib inflateInit failed: " << (int)status << std::endl;
		exit(1);
	}

	// decompress
	status = inflate(&zs, Z_FINISH);
	if(status != Z_STREAM_END){
		inflateEnd(&zs);
		std::cerr << Helpers::timestamp("ERROR","BGZF") << "Zlib inflateEnd failed: " << (int)status << std::endl;
		exit(1);
	}

	// finalize
	status = inflateEnd(&zs);
	if(status != Z_OK){
		inflateEnd(&zs);
		std::cerr << Helpers::timestamp("ERROR","BGZF") << "Zlib inflateFinalize failed: " << (int)status << std::endl;
		exit(1);
	}

	//if(zs.total_out == 0)
	//	std::cerr << Helpers::timestamp("LOG", "BGZF") << "Detected empty BGZF block" << std::endl;

	output.pointer += zs.total_out;

	return(true);
}

} /* namespace IO */
} /* namespace Tomahawk */
