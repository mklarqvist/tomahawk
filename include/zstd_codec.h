#ifndef ALGORITHM_COMPRESSION_ZSTD_CODEC_H_
#define ALGORITHM_COMPRESSION_ZSTD_CODEC_H_

#include <cstdint>

#include "buffer.h"

namespace tomahawk {

class ZSTDCodec{
public:
	class ZSTDCodecInternal;

	ZSTDCodec();
	~ZSTDCodec();

	ZSTDCodec(const ZSTDCodec& other) = delete;
	ZSTDCodec& operator=(const ZSTDCodec& other) = delete;

	bool InitStreamCompress(const int compression_level);
	bool StopStreamCompress();
	bool InitStreamDecompress();
	size_t StreamCompress(const twk_buffer_t& src, twk_buffer_t& dst, std::ostream& out, const uint32_t block_size = 512000);
	bool StreamDecompress(twk_buffer_t& src, twk_buffer_t& dst);

	bool Compress(const twk_buffer_t& src, twk_buffer_t& dst, const int compression_level);
	bool Decompress(const twk_buffer_t& src, twk_buffer_t& dst);

	std::ostream& WriteOutbuf(std::ostream& ostream);
	size_t GetOutputSize() const;

public:
	ZSTDCodecInternal* mImpl;
};

}



#endif /* ALGORITHM_COMPRESSION_ZSTD_CODEC_H_ */
