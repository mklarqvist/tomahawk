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
#ifndef TWK_READER_H_
#define TWK_READER_H_

#include <fstream>

#include "core.h"
#include "header.h"
#include "index.h"
#include "zstd_codec.h"

namespace tomahawk {

/**<
 * Block iterator for tomahawk blocks (twk). Allows the iteration of blocks
 * either as raw compressed data blocks or pre-processed uncompressed blocks.
 */
class twk1_blk_iterator {
public:
	twk1_blk_iterator() : stream(nullptr){}
	~twk1_blk_iterator(){ }

	bool NextBlockRaw();
	bool NextBlock();
	inline const twk1_block_t& GetBlock(void) const{ return(this->blk); }

public:
	ZSTDCodec zcodec; // support codec
	twk_buffer_t buf; // support buffer
	twk_oblock_t oblk;// block wrapper
	twk1_block_t blk; // block
	std::istream* stream; // stream pointer
};

/**<
 * Reader of twk files.
 */
class twk_reader {
public:
	twk_reader() : buf(nullptr), stream(nullptr){}
	~twk_reader(){ delete stream; }

	/**<
	 * Open a target twk file. File header, index, and footer will be read
	 * and parsed as part of the opening procedure. If these pass without errors
	 * then return TRUE. Returns FALSE otherwise.
	 * @param file Input target twk file.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool Open(std::string file);

public:
	std::streambuf* buf;
	std::istream*   stream;
	std::ifstream   fstream;
	VcfHeader   hdr;
	Index index;
};

}



#endif /* TWK_READER_H_ */
