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

	bool Compress(const twk_buffer_t& src, twk_buffer_t& dst, const int compression_level);
	bool Decompress(const twk_buffer_t& src, twk_buffer_t& dst);

private:
	ZSTD_CCtx* compression_context_; // recycle contexts
	ZSTD_DCtx* decompression_context_; // recycle contexts
};

}



#endif /* ALGORITHM_COMPRESSION_ZSTD_CODEC_H_ */
