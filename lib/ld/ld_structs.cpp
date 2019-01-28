#include "twk_reader.h"
#include "ld_structs.h"

namespace tomahawk {

/****************************
*  twk1_ldd_blk
****************************/
twk1_ldd_blk::twk1_ldd_blk() : owns_block(false), unphased(true), n_rec(0), m_vec(0), m_list(0), m_bitmap(0), blk(nullptr), vec(nullptr), list(nullptr), bitmap(nullptr){}
twk1_ldd_blk::twk1_ldd_blk(twk1_blk_iterator& it, const uint32_t n_samples) :
	owns_block(false), unphased(true),
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
	unphased = other.unphased;
	owns_block = false; n_rec = 0; // do not change m
	blk = other.blk;
	if(blk != nullptr){
		n_rec = other.blk->n;
	}
	return(*this);
}

twk1_ldd_blk& twk1_ldd_blk::operator=(twk1_ldd_blk&& other){
	if(owns_block) delete blk;
	unphased = other.unphased;
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
	unphased = other.unphased;
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
		// If we have alrady unpacked the genotypes into a bitvector
		// in the step above then this list object do not own this
		// object.
		if(unpack & TWK_LDD_VEC){
			for(int i = 0; i < blk->n; ++i){
				if(blk->rcds[i].an) continue; // do not construct if missing data
				list[i].own = false;
				list[i].n   = vec[i].n;
				list[i].bv  = vec[i].data;
				list[i].Build(blk->rcds[i], n_samples, resizeable, false);
			}
		} else {
			for(int i = 0; i < blk->n; ++i){
				if(blk->rcds[i].an) continue; // do not construct if missing data
				list[i].Build(blk->rcds[i], n_samples, resizeable, false);
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

}
