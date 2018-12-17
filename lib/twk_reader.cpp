#include "twk_reader.h"

namespace tomahawk {

/****************************
*  twk1_blk_iterator
****************************/
bool twk1_blk_iterator::NextBlockRaw(){
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

bool twk1_blk_iterator::NextBlock(){
	if(this->NextBlockRaw() == false)
		return false;

	// Decompress data
	zcodec.Decompress(oblk.bytes, buf);
	buf >> blk;
	buf.reset();

	return true;
}

/****************************
*  twk1_ldd_blk
****************************/
twk1_ldd_blk::twk1_ldd_blk() : owns_block(false), n_rec(0), m_vec(0), m_list(0), m_bitmap(0), blk(nullptr), vec(nullptr), list(nullptr), bitmap(nullptr){}
twk1_ldd_blk::twk1_ldd_blk(twk1_blk_iterator& it, const uint32_t n_samples) :
	owns_block(false),
	n_rec(it.blk.n),
	m_vec(0), m_list(0), m_bitmap(0),
	blk(&it.blk),
	vec(new twk_igt_vec[it.blk.n]),
	list(new twk_igt_list[it.blk.n]),
	bitmap(nullptr)
{
	for(int i = 0; i < it.blk.n; ++i){
		vec[i].Build(it.blk.rcds[i], n_samples);
		list[i].Build(it.blk.rcds[i], n_samples);
	}
}

twk1_ldd_blk::~twk1_ldd_blk(){
	if(owns_block) delete blk;
	delete[] vec;
	delete[] list;
	delete[] bitmap;
}

twk1_ldd_blk& twk1_ldd_blk::operator=(const twk1_ldd_blk& other){
	if(owns_block) delete blk;
	owns_block = false; n_rec = 0; // do not change m
	blk = other.blk;
	if(blk != nullptr){
		n_rec = other.blk->n;
	}
	return(*this);
}

twk1_ldd_blk& twk1_ldd_blk::operator=(twk1_ldd_blk&& other){
	if(owns_block) delete blk;
	blk = nullptr;
	std::swap(blk, other.blk);
	owns_block = other.owns_block;
	n_rec = other.n_rec; // do not change m

	delete[] vec; vec = nullptr;
	delete[] list; list = nullptr;
	delete[] bitmap; bitmap = nullptr;
	std::swap(vec, other.vec);
	std::swap(list, other.list);
	std::swap(bitmap, other.bitmap);

	if(blk != nullptr){
		n_rec = blk->n;
	}

	return(*this);
}

void twk1_ldd_blk::SetPreloaded(const twk1_ldd_blk& other){
	if(owns_block) delete blk;
	owns_block = false; n_rec = 0; // do not change m
	blk = other.blk;
	if(blk != nullptr){
		n_rec = other.blk->n;
	}
	this->list = other.list;
	this->vec = other.vec;
	this->m_list = other.m_list;
	this->m_vec = other.m_vec;
	this->bitmap = other.bitmap;
	this->m_bitmap = other.m_bitmap;
}

void twk1_ldd_blk::Set(twk1_blk_iterator& it, const uint32_t n_samples){
	if(owns_block) delete blk;
	delete[] vec;
	delete[] list;

	owns_block = false;
	n_rec = it.blk.n;
	blk = &it.blk;
}

void twk1_ldd_blk::SetOwn(twk1_blk_iterator& it, const uint32_t n_samples){
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

void twk1_ldd_blk::SetOwn(twk1_block_t& it, const uint32_t n_samples){
	if(owns_block) delete blk;
	delete[] vec; vec = nullptr;
	delete[] list; list = nullptr;

	owns_block = false;
	blk = &it;
	// Decompress data
	//it.zcodec.Decompress(it.oblk.bytes, it.buf);
	//it.buf >> *blk;
	//it.buf.reset();
	n_rec = blk->n;
}


void twk1_ldd_blk::Clear(){
	delete[] vec;
	delete[] list;
}

void twk1_ldd_blk::Inflate(const uint32_t n_samples,
                           const uint8_t unpack,
                           const bool resizeable)
{
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

	if(unpack & TWK_LDD_BITMAP){
		if(blk->n > m_bitmap){
			bitmap = new bitmap_type[blk->n];
			m_bitmap = blk->n;
		}
	}

	if(unpack & TWK_LDD_VEC){
		for(int i = 0; i < blk->n; ++i)
			vec[i].Build(blk->rcds[i], n_samples);
	}

	if(unpack & TWK_LDD_LIST){
		if(unpack & TWK_LDD_VEC){
			//std::cerr << "not owner" << std::endl;
			for(int i = 0; i < blk->n; ++i){
				if(blk->rcds[i].an) continue; // do not construct if missing data
				list[i].own = false; list[i].n = vec[i].n;
				list[i].bv = vec[i].data;
				list[i].Build(blk->rcds[i], n_samples, resizeable);
			}
		} else {
			for(int i = 0; i < blk->n; ++i){
				if(blk->rcds[i].an) continue; // do not construct if missing data
				list[i].Build(blk->rcds[i], n_samples, resizeable);
			}
		}

	}

	if(unpack & TWK_LDD_BITMAP){
		for(int i = 0; i < blk->n; ++i){
			if(blk->rcds[i].an) continue; // do not construct if missing data
			bitmap[i].reset(); // does not release memory used
			uint32_t cumpos = 0;
			// iterate over run-length encoded entries
			uint32_t k = 0;
			for(int j = 0; j < blk->rcds[i].gt->n; ++j){
				const uint32_t len  = blk->rcds[i].gt->GetLength(j);
				const uint8_t  refA = blk->rcds[i].gt->GetRefA(j);
				const uint8_t  refB = blk->rcds[i].gt->GetRefB(j);

				if(refA == 0 && refB == 0){
					cumpos += 2*len;
					continue;
				}

				for(int k = 0; k < 2*len; k+=2){
					if(refA != 0){ bitmap[i].set(cumpos + k + 0); }
					if(refB != 0){ bitmap[i].set(cumpos + k + 1); }
				}
				cumpos += 2*len;
			}
			//	std::cerr << bitmap[i] << std::endl;
		}
	}
}

/****************************
*  twk_reader
****************************/
bool twk_reader::Open(std::string file){
	fstream.open(file, std::ios::in|std::ios::binary|std::ios::ate);
	if(!fstream.good()){
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to open \"" << file << "\"!" << std::endl;
		return false;
	}
	buf = fstream.rdbuf();
	stream = new std::istream(buf);

	//stream->open(file, std::ios::in|std::ios::binary|std::ios::ate);
	uint64_t filesize = stream->tellg();
	stream->seekg(0);

	// read magic
	char magic[TOMAHAWK_MAGIC_HEADER_LENGTH];
	stream->read(magic, TOMAHAWK_MAGIC_HEADER_LENGTH);
	if(strncmp(magic, TOMAHAWK_MAGIC_HEADER.data(), TOMAHAWK_MAGIC_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to read MAGIC!" << std::endl;
		return false;
	}

	uint64_t buf_size = 0, obuf_size = 0;
	stream->read(reinterpret_cast<char*>(&buf_size), sizeof(uint64_t));
	stream->read(reinterpret_cast<char*>(&obuf_size),sizeof(uint64_t));
	twk_buffer_t obuf(obuf_size);
	twk_buffer_t buf(buf_size);
	stream->read(obuf.data(),obuf_size);
	obuf.n_chars_ = obuf_size;
	//std::cerr << "header=" << buf_size << "," << obuf_size << "/" << buf.capacity() << "/" << obuf.capacity() << std::endl;

	ZSTDCodec zcodec;
	if(zcodec.Decompress(obuf, buf) == false){
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to decompress header!" << std::endl;
		return false;
	}
	//std::cerr << "bufs=" << buf.size() << "==" << buf_size << std::endl;
	assert(buf.size() == buf_size);
	buf >> hdr;
	buf.reset(); obuf.reset();
	//std::cerr << "done hdr" << std::endl;

	uint64_t data_start = stream->tellg();

	// seek to end-of-file
	stream->seekg(filesize - TOMAHAWK_FILE_EOF_LENGTH - sizeof(uint64_t));
	uint64_t offset_start_index = 0;
	stream->read(reinterpret_cast<char*>(&offset_start_index), sizeof(uint64_t));
	//std::cerr << "seek offset=" << offset_start_index << "/" << filesize << std::endl;
	stream->seekg(offset_start_index);
	if(stream->good() == false){
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to seek in file!" << std::endl;
		return false;
	}
	//std::cerr << "seek good=" << stream->tellg() << "/" << filesize << std::endl;

	uint8_t marker = 0;
	stream->read(reinterpret_cast<char*>(&marker),   sizeof(uint8_t));
	stream->read(reinterpret_cast<char*>(&buf_size), sizeof(uint64_t));
	stream->read(reinterpret_cast<char*>(&obuf_size),sizeof(uint64_t));
	obuf.resize(obuf_size), buf.resize(buf_size);
	//std::cerr << "before read=" << obuf_size << std::endl;
	stream->read(obuf.data(),obuf_size);
	obuf.n_chars_ = obuf_size;
	//std::cerr << "header=" << buf_size << "," << obuf_size << "/" << buf.capacity() << "/" << obuf.capacity() << std::endl;


	if(zcodec.Decompress(obuf, buf) == false){
		std::cerr << utility::timestamp("ERROR","TWK") << "Failed to decompress index!" << std::endl;
		return false;
	}
	buf >> index;

	// Seek back to start of data.
	stream->seekg(data_start);

	return true;
}

}
