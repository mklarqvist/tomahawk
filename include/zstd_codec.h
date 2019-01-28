/*
Copyright (C) 2016-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TWK_ZSTD_CODEC_H_
#define TWK_ZSTD_CODEC_H_

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

#endif /* TWK_ZSTD_CODEC_H_ */
