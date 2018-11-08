#include "zstd_codec.h"

namespace tomahawk {

ZSTDCodec::ZSTDCodec() :
	compression_context_(ZSTD_createCCtx()),
	decompression_context_(ZSTD_createDCtx()),
	cstream_context(ZSTD_createCStream()),
	dstream_context(ZSTD_createDStream())
{

}

ZSTDCodec::~ZSTDCodec(){
	ZSTD_freeCCtx(compression_context_);
	ZSTD_freeDCtx(decompression_context_);
	ZSTD_freeCStream(cstream_context);
	ZSTD_freeDStream(dstream_context);
}

bool ZSTDCodec::InitStreamCompress(const int compression_level){
	ZSTD_initCStream(cstream_context, compression_level);
	return true;
}
bool ZSTDCodec::StopStreamCompress(){
	assert(ZSTD_flushStream(cstream_context, &outbuf) == 0);
	assert(ZSTD_endStream(cstream_context, &outbuf) == 0);
	return true;
}

bool ZSTDCodec::InitStreamDecompress(){
	ZSTD_initDStream(dstream_context);
	return true;
}

size_t ZSTDCodec::StreamCompress(const twk_buffer_t& src,
		twk_buffer_t& dst,
		std::ostream& out,
		const uint32_t block_size)
	{
	size_t left = src.size();
	dst.resize(block_size + 65536);

	outbuf.dst  = dst.data();
	outbuf.size = dst.capacity();
	outbuf.pos  = 0;

	const char* src_pos = src.data();
	size_t n_c = 0;

	while (left) {
		inbuf.pos  = 0;
		inbuf.size = (left > block_size ? block_size : left);
		inbuf.src  = src_pos;

		const size_t ret = ZSTD_compressStream(cstream_context, &outbuf, &inbuf);
		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			return(-1);
		}

		ZSTD_flushStream(cstream_context, &outbuf);
		out.write((const char*)outbuf.dst, outbuf.pos);
		//std::cerr << "wrote=" << inbuf.pos << "->" << outbuf.pos << " = " << (float)inbuf.pos/outbuf.pos << "-fold" << std::endl;
		n_c += outbuf.pos;
		outbuf.pos = 0;

		src_pos += inbuf.size; // move source pos up
		left    -= inbuf.size; // decrement remainder
	}
	return(n_c);
}

bool ZSTDCodec::StreamDecompress(twk_buffer_t& src, twk_buffer_t& dst){
	dst.resize(src.size()	);

	outbuf.dst  = dst.data();
	outbuf.size = dst.capacity();
	outbuf.pos  = dst.n_chars_;

	inbuf.pos  = src.iterator_position_;
	inbuf.size = src.size();
	inbuf.src  = src.data();

	//const char* src_pos = src.data();

	while(true){
		const size_t ret = ZSTD_decompressStream(dstream_context, &outbuf, &inbuf);
		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			return(false);
		}
		dst.n_chars_ = outbuf.pos;
		if(inbuf.pos == inbuf.size) break;
		dst.resize(dst.capacity() * 2);
		outbuf.size = dst.capacity();
		outbuf.dst = dst.data();
	}

	src.iterator_position_ = inbuf.size;
	return(inbuf.pos == inbuf.size);

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
