#ifndef TWK_READER_H_
#define TWK_READER_H_

#include <fstream>

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

#define TWK_LDD_NONE 0
#define TWK_LDD_VEC  1
#define TWK_LDD_LIST 2
#define TWK_LDD_ALL  ((TWK_LDD_VEC) | (TWK_LDD_LIST))

struct twk1_ldd_blk {
	twk1_ldd_blk() : owns_block(false), n_rec(0), m_vec(0), m_list(0), blk(nullptr), vec(nullptr), list(nullptr){}
	twk1_ldd_blk(twk1_blk_iterator& it, const uint32_t n_samples) :
		owns_block(false),
		n_rec(it.blk.n),
		m_vec(0), m_list(0),
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

	twk1_ldd_blk& operator=(const twk1_ldd_blk& other){
		if(owns_block) delete blk;
		owns_block = false; n_rec = 0; // do not change m
		blk = other.blk;
		if(blk != nullptr){
			n_rec = other.blk->n;
		}
		return(*this);
	}

	void Set(twk1_blk_iterator& it, const uint32_t n_samples){
		if(owns_block) delete blk;
		delete[] vec;
		delete[] list;

		owns_block = false;
		n_rec = it.blk.n;
		blk = &it.blk;
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
	}

	void Clear(){
		delete[] vec;
		delete[] list;
	}

	void Inflate(const uint32_t n_samples, const uint8_t unpack = TWK_LDD_ALL){
		if(unpack & TWK_LDD_VEC){
			if(blk->n > m_vec){
				delete[] vec;
				vec  = new twk_igt_vec[blk->n];
				m_vec = blk->n;
			}
		}
		if(unpack & TWK_LDD_LIST){
			if(blk->n > m_list){
				list = new twk_igt_list[blk->n];
				m_list = blk->n;
			}
		}

		if(unpack & TWK_LDD_VEC){
			for(int i = 0; i < blk->n; ++i)
				vec[i].Build(blk->rcds[i], n_samples);
		}

		if(unpack & TWK_LDD_LIST){
			for(int i = 0; i < blk->n; ++i)
				list[i].Build(blk->rcds[i], n_samples);
		}
	}

	inline void operator=(twk1_block_t* block){ this->blk = block; }
	inline void operator=(twk1_block_t& block){ this->blk = &block; }
	inline const twk1_t& operator[](const uint32_t p) const{ return(blk->rcds[p]); }

	bool owns_block;
	uint32_t n_rec, m_vec, m_list;
	twk1_block_t* blk; // data block
	twk_igt_vec* vec; // vectorized (bitvector)
	twk_igt_list* list; // list
};

/**<
 * Basic reader
 */
class twk_reader {
public:
	twk_reader() : rstream(nullptr){}
	~twk_reader(){
		delete rstream;
	}

	bool Open(void);
	bool Open(std::string file){
		rstream = new std::ifstream;
		std::ifstream* stream = reinterpret_cast<std::ifstream*>(rstream);
		stream->open(file, std::ios::in|std::ios::binary|std::ios::ate);
		if(!stream->good()){
			std::cerr << "failed to open: " << file << std::endl;
			return false;
		}
		uint64_t filesize = stream->tellg();
		stream->seekg(0);

		// read magic
		char magic[TOMAHAWK_MAGIC_HEADER_LENGTH];
		stream->read(magic, TOMAHAWK_MAGIC_HEADER_LENGTH);
		if(strncmp(magic, TOMAHAWK_MAGIC_HEADER.data(), TOMAHAWK_MAGIC_HEADER_LENGTH) != 0){
			std::cerr << "failed to read maagic" << std::endl;
			return false;
		}

		uint64_t buf_size = 0, obuf_size = 0;
		stream->read(reinterpret_cast<char*>(&buf_size), sizeof(uint64_t));
		stream->read(reinterpret_cast<char*>(&obuf_size),sizeof(uint64_t));
		twk_buffer_t obuf(obuf_size);
		twk_buffer_t buf(buf_size);
		stream->read(obuf.data(),obuf_size);
		obuf.n_chars_ = obuf_size;
		std::cerr << "header=" << buf_size << "," << obuf_size << "/" << buf.capacity() << "/" << obuf.capacity() << std::endl;

		ZSTDCodec zcodec;
		if(zcodec.Decompress(obuf, buf) == false){
			std::cerr << "failed to decompress header" << std::endl;
			return false;
		}
		std::cerr << "bufs=" << buf.size() << "==" << buf_size << std::endl;
		assert(buf.size() == buf_size);
		buf >> hdr;
		buf.reset(); obuf.reset();
		std::cerr << "done hdr" << std::endl;

		uint64_t data_start = stream->tellg();

		// seek to end-of-file
		stream->seekg(filesize - TOMAHAWK_FILE_EOF_LENGTH - sizeof(uint64_t));
		uint64_t offset_start_index = 0;
		stream->read(reinterpret_cast<char*>(&offset_start_index), sizeof(uint64_t));
		std::cerr << "seek offset=" << offset_start_index << "/" << filesize << std::endl;
		stream->seekg(offset_start_index);
		if(stream->good() == false){
			std::cerr << "failed seek" << std::endl;
		}
		std::cerr << "seek good=" << stream->tellg() << "/" << filesize << std::endl;

		uint8_t marker = 0;
		stream->read(reinterpret_cast<char*>(&marker),   sizeof(uint8_t));
		stream->read(reinterpret_cast<char*>(&buf_size), sizeof(uint64_t));
		stream->read(reinterpret_cast<char*>(&obuf_size),sizeof(uint64_t));
		obuf.resize(obuf_size), buf.resize(buf_size);
		std::cerr << "before read=" << obuf_size << std::endl;
		stream->read(obuf.data(),obuf_size);
		obuf.n_chars_ = obuf_size;
		std::cerr << "header=" << buf_size << "," << obuf_size << "/" << buf.capacity() << "/" << obuf.capacity() << std::endl;


		if(zcodec.Decompress(obuf, buf) == false){
			std::cerr << "failed to decompress" << std::endl;
			return false;
		}
		buf >> index;

		stream->seekg(data_start);

		return true;
	}

public:
	std::istream* rstream;
	io::VcfHeader hdr;
	Index index;
};

}



#endif /* TWK_READER_H_ */
