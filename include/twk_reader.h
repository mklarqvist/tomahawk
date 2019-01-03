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
