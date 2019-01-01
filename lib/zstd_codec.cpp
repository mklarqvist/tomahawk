#include "zstd.h"
#include "zstd_errors.h"
#include "zstd_codec.h"

namespace tomahawk {

class ZSTDCodec::ZSTDCodecInternal {
private:
	typedef        ZSTDCodec   self_type;
	typedef struct ZSTD_CCtx_s ZSTD_CCtx;
	typedef struct ZSTD_DCtx_s ZSTD_DCtx;

public:
	ZSTDCodecInternal();
	~ZSTDCodecInternal();

	std::ostream& WriteOutbuf(std::ostream& ostream){
		ostream.write((const char*)outbuf.dst, outbuf.pos);
		return(ostream);
	}

	size_t GetOutputSize() const { return(outbuf.pos); }

	ZSTD_inBuffer_s inbuf;
	ZSTD_outBuffer_s outbuf;
	ZSTD_CCtx* compression_context_; // recycle contexts
	ZSTD_DCtx* decompression_context_; // recycle contexts
	ZSTD_CStream* cstream_context;
	ZSTD_DStream* dstream_context;
};

ZSTDCodec::ZSTDCodec() : mImpl(new ZSTDCodecInternal){}
ZSTDCodec::~ZSTDCodec(){ delete mImpl; }

ZSTDCodec::ZSTDCodecInternal::ZSTDCodecInternal() :
	compression_context_(ZSTD_createCCtx()),
	decompression_context_(ZSTD_createDCtx()),
	cstream_context(ZSTD_createCStream()),
	dstream_context(ZSTD_createDStream())
{

}

ZSTDCodec::ZSTDCodecInternal::~ZSTDCodecInternal(){
	ZSTD_freeCCtx(compression_context_);
	ZSTD_freeDCtx(decompression_context_);
	ZSTD_freeCStream(cstream_context);
	ZSTD_freeDStream(dstream_context);
}

bool ZSTDCodec::InitStreamCompress(const int compression_level){
	ZSTD_initCStream(mImpl->cstream_context, compression_level);
	return true;
}
bool ZSTDCodec::StopStreamCompress(){
	assert(ZSTD_flushStream(mImpl->cstream_context, &mImpl->outbuf) == 0);
	assert(ZSTD_endStream(mImpl->cstream_context, &mImpl->outbuf) == 0);
	return true;
}

bool ZSTDCodec::InitStreamDecompress(){
	ZSTD_initDStream(mImpl->dstream_context);
	return true;
}

size_t ZSTDCodec::StreamCompress(const twk_buffer_t& src,
		twk_buffer_t& dst,
		std::ostream& out,
		const uint32_t block_size)
	{
	size_t left = src.size();
	dst.resize(block_size + 65536);

	mImpl->outbuf.dst  = dst.data();
	mImpl->outbuf.size = dst.capacity();
	mImpl->outbuf.pos  = 0;

	const char* src_pos = src.data();
	size_t n_c = 0;

	while (left) {
		mImpl->inbuf.pos  = 0;
		mImpl->inbuf.size = (left > block_size ? block_size : left);
		mImpl->inbuf.src  = src_pos;

		const size_t ret = ZSTD_compressStream(mImpl->cstream_context, &mImpl->outbuf, &mImpl->inbuf);
		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			return(-1);
		}

		ZSTD_flushStream(mImpl->cstream_context, &mImpl->outbuf);
		out.write((const char*)mImpl->outbuf.dst, mImpl->outbuf.pos);
		//std::cerr << "wrote=" << inbuf.pos << "->" << outbuf.pos << " = " << (float)inbuf.pos/outbuf.pos << "-fold" << std::endl;
		n_c += mImpl->outbuf.pos;
		mImpl->outbuf.pos = 0;

		src_pos += mImpl->inbuf.size; // move source pos up
		left    -= mImpl->inbuf.size; // decrement remainder
	}
	return(n_c);
}

bool ZSTDCodec::StreamDecompress(twk_buffer_t& src, twk_buffer_t& dst){
	dst.resize(src.size()	);

	mImpl->outbuf.dst  = dst.data();
	mImpl->outbuf.size = dst.capacity();
	mImpl->outbuf.pos  = dst.n_chars_;

	mImpl->inbuf.pos  = src.iterator_position_;
	mImpl->inbuf.size = src.size();
	mImpl->inbuf.src  = src.data();

	//const char* src_pos = src.data();

	while(true){
		const size_t ret = ZSTD_decompressStream(mImpl->dstream_context, &mImpl->outbuf, &mImpl->inbuf);
		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			return(false);
		}
		dst.n_chars_ = mImpl->outbuf.pos;
		if(mImpl->inbuf.pos == mImpl->inbuf.size) break;
		dst.resize(dst.capacity() * 2);
		mImpl->outbuf.size = dst.capacity();
		mImpl->outbuf.dst = dst.data();
	}

	src.iterator_position_ = mImpl->inbuf.size;
	return(mImpl->inbuf.pos == mImpl->inbuf.size);

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

std::ostream& ZSTDCodec::WriteOutbuf(std::ostream& ostream){
	return(mImpl->WriteOutbuf(ostream));
}

size_t ZSTDCodec::GetOutputSize() const{
	return(mImpl->GetOutputSize());
}

}
