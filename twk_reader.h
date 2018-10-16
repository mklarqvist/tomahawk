#ifndef TWK_READER_H_
#define TWK_READER_H_

#include "vcf_utils.h"
#include "index.h"
#include "zstd_codec.h"

namespace tomahawk {

/**<
 * Tomahawk block iterator.
 */
class twk1_blk_iterator {
public:
	twk1_blk_iterator() : stream(nullptr){}
	~twk1_blk_iterator(){ }

	bool NextBlockRaw(){
		if(stream->good() == false){
			std::cerr << "stream died" << std::endl;
			return false;
		}

		uint8_t marker = 0;
		DeserializePrimitive(marker, *stream);
		if(marker == 0){
			//std::cerr << "0 marker found. stopping" << std::endl;
			return false;
		}
		assert(marker == 1);

		*stream >> oblk;
		if(stream->good() == false){
			std::cerr << "stream died" << std::endl;
			return false;
		}

		assert(oblk.bytes.size() == oblk.nc);
		buf.resize(oblk.n);

		return true;
	}

	bool NextBlock(){
		if(this->NextBlockRaw() == false)
			return false;

		// Decompress data
		zcodec.Decompress(oblk.bytes, buf);
		buf >> blk;
		buf.reset();

		return true;
	}

	inline const twk1_block_t& GetBlock(void) const{ return(this->blk); }

public:
	ZSTDCodec zcodec; // support codec
	twk_buffer_t buf; // support buffer
	twk_oblock_t oblk;// block wrapper
	twk1_block_t blk; // block
	std::istream* stream; // stream pointer
};

struct twk1_ldd_blk {
	twk1_ldd_blk() : owns_block(0), n_rec(0), blk(nullptr), vec(nullptr), list(nullptr){}
	twk1_ldd_blk(twk1_blk_iterator& it, const uint32_t n_samples) :
		owns_block(false),
		n_rec(it.blk.n),
		blk(&it.blk),
		vec(new twk_igt_vec[it.blk.n]),
		list(new twk_igt_list[it.blk.n])
	{
		for(int i = 0; i < it.blk.n; ++i){
			vec[i].Build(it.blk.rcds[i], n_samples);
			list[i].Build(it.blk.rcds[i], n_samples);
		}
	}

	~twk1_ldd_blk(){
		if(owns_block) delete blk;
		delete[] vec;
		delete[] list;
	}

	void Set(twk1_blk_iterator& it, const uint32_t n_samples){
		if(owns_block) delete blk;
		delete[] vec;
		delete[] list;

		owns_block = false;
		n_rec = it.blk.n;
		blk = &it.blk;
		vec = new twk_igt_vec[it.blk.n];
		list = new twk_igt_list[it.blk.n];

		for(int i = 0; i < it.blk.n; ++i){
			vec[i].Build(it.blk.rcds[i], n_samples);
			list[i].Build(it.blk.rcds[i], n_samples);
		}
	}

	void SetOwn(twk1_blk_iterator& it, const uint32_t n_samples){
		if(owns_block) delete blk;
		delete[] vec; vec = nullptr;
		delete[] list; list = nullptr;

		owns_block = true;
		blk = new twk1_block_t;
		// Decompress data
		it.zcodec.Decompress(it.oblk.bytes, it.buf);
		it.buf >> *blk;
		it.buf.reset();
		n_rec = blk->n;

		/*
		vec  = new twk_igt_vec[blk->n];
		list = new twk_igt_list[blk->n];

		for(int i = 0; i < blk->n; ++i){
			vec[i].Build(blk->rcds[i], n_samples);
			list[i].Build(blk->rcds[i], n_samples);
		}
		*/
	}

	void Clear(){
		delete[] vec;
		delete[] list;
	}

	void Inflate(const uint32_t n_samples){
		vec  = new twk_igt_vec[blk->n];
		list = new twk_igt_list[blk->n];

		for(int i = 0; i < blk->n; ++i){
			vec[i].Build(blk->rcds[i], n_samples);
			list[i].Build(blk->rcds[i], n_samples);
		}
	}

	inline void operator=(twk1_block_t* block){ this->blk = block; }
	inline void operator=(twk1_block_t& block){ this->blk = &block; }
	inline const twk1_t& operator[](const uint32_t p) const{ return(blk->rcds[p]); }

	bool owns_block;
	uint32_t n_rec;
	twk1_block_t* blk; // data block
	twk_igt_vec* vec; // vectorized (bitvector)
	twk_igt_list* list; // list
};

/**<
 * Basic reader
 */
class twk_reader {
public:

	bool Open(void);
	bool Open(std::string file);

public:
	io::VcfHeader hdr;
	Index index;
};

}



#endif /* TWK_READER_H_ */
