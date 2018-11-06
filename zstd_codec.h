#ifndef ALGORITHM_COMPRESSION_ZSTD_CODEC_H_
#define ALGORITHM_COMPRESSION_ZSTD_CODEC_H_

#include "zstd.h"
#include "zstd_errors.h"
#include <cstdint>

#include "buffer.h"

namespace tomahawk {

class ZSTDCodec{
private:
	typedef        ZSTDCodec   self_type;
	typedef struct ZSTD_CCtx_s ZSTD_CCtx;
	typedef struct ZSTD_DCtx_s ZSTD_DCtx;

public:
	ZSTDCodec();
	~ZSTDCodec();

	bool InitStreamCompress(const int compression_level);
	bool StopStreamCompress();
	bool InitStreamDecompress();
	size_t StreamCompress(const twk_buffer_t& src, twk_buffer_t& dst, std::ostream& out, const uint32_t block_size = 512000);
	bool StreamDecompress(twk_buffer_t& src, twk_buffer_t& dst);

	bool Compress(const twk_buffer_t& src, twk_buffer_t& dst, const int compression_level);
	bool Decompress(const twk_buffer_t& src, twk_buffer_t& dst);

public:
	ZSTD_inBuffer_s inbuf;
	ZSTD_outBuffer_s outbuf;
	ZSTD_CCtx* compression_context_; // recycle contexts
	ZSTD_DCtx* decompression_context_; // recycle contexts
	ZSTD_CStream* cstream_context;
	ZSTD_DStream* dstream_context;
};

}



#endif /* ALGORITHM_COMPRESSION_ZSTD_CODEC_H_ */
