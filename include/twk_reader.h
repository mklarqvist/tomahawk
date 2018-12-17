#ifndef TWK_READER_H_
#define TWK_READER_H_

#include <fstream>

#include "core.h"
#include "vcf_utils.h"
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
 * Expanded twk block for use in linkage-disequilibrium calculations. Expands
 * data into:
 *    1) bit-index lists;
 *    2) uncompressed bitvectors;
 *    3) compressed run-length encoded genotypes;
 *    4) compressed hybrid bitmaps (EWAH);
 * What data to pre-process is determined by passing the appropriate `TWK_LDD_`-prefixed
 * macro-definied value: TWK_LDD_NONE, TWK_LDD_VEC, TWK_LDD_LIST, and TWK_LDD_ALL.
 */
struct twk1_ldd_blk {
	typedef EWAHBoolArray<uint64_t> bitmap_type;

	twk1_ldd_blk();
	twk1_ldd_blk(twk1_blk_iterator& it, const uint32_t n_samples);
	~twk1_ldd_blk();
	twk1_ldd_blk& operator=(const twk1_ldd_blk& other);
	twk1_ldd_blk& operator=(twk1_ldd_blk&& other);

	void SetPreloaded(const twk1_ldd_blk& other);
	void Set(twk1_blk_iterator& it, const uint32_t n_samples);
	void SetOwn(twk1_blk_iterator& it, const uint32_t n_samples);
	void SetOwn(twk1_block_t& it, const uint32_t n_samples);
	void Clear();

	void Inflate(const uint32_t n_samples,
	             const uint8_t unpack = TWK_LDD_ALL,
	             const bool resizeable = false);

	inline void operator=(twk1_block_t* block){ this->blk = block; }
	inline void operator=(twk1_block_t& block){ this->blk = &block; }
	inline const twk1_t& operator[](const uint32_t p) const{ return(blk->rcds[p]); }

public:
	bool owns_block;
	uint32_t n_rec, m_vec, m_list, m_bitmap; // number records, memory allocated for vectors, memory allocated for lists
	twk1_block_t* blk; // data block
	twk_igt_vec*  vec; // vectorized (bitvector)
	twk_igt_list* list; // list
	bitmap_type* bitmap; // bitmap
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
	io::VcfHeader   hdr;
	Index index;
};

}



#endif /* TWK_READER_H_ */
