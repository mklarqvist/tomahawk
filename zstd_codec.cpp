#include "zstd_codec.h"

namespace tomahawk {

ZSTDCodec::ZSTDCodec() :
	compression_context_(ZSTD_createCCtx()),
	decompression_context_(ZSTD_createDCtx())
{

}

ZSTDCodec::~ZSTDCodec(){
	ZSTD_freeCCtx(this->compression_context_);
	ZSTD_freeDCtx(this->decompression_context_);
}

bool ZSTDCodec::Compress(const twk_buffer_t& src, twk_buffer_t& dst, const int compression_level){
	dst.reset();
	dst.resize(src.size() + 65536);
	const size_t ret = ZSTD_compress( dst.data(), dst.capacity(),
	                                  src.data(), src.size(),
	                                  compression_level);

	//std::cerr << utility::timestamp("LOG","COMPRESSION") << "Input: " << src.size() << " and output: " << ret << " -> " << (float)src.size()/ret << "-fold"  << std::endl;

	if(ZSTD_isError(ret)){
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		return(false);
	}
	dst.n_chars_ = ret;

	return true;
}

bool ZSTDCodec::Decompress(const twk_buffer_t& src, twk_buffer_t& dst){
	const size_t ret = ZSTD_decompress( dst.data(), dst.capacity(),
	                                    src.data(), src.size());

	//std::cerr << utility::timestamp("LOG","COMPRESSION") << "Input: " << src.size() << " and output: " << ret << " -> " << (float)ret/src.size() << "-fold"  << std::endl;

	if(ZSTD_isError(ret)){
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		return(false);
	}

	dst.n_chars_ = ret;

	return true;
}

}
