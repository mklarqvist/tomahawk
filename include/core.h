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
#ifndef TWK_CORE_H_
#define TWK_CORE_H_

#include <cstdint>

#include "tomahawk.h"
#include "buffer.h"
#include "third_party/ewah.h"

namespace tomahawk {

/****************************
*  Core support
****************************/
const uint8_t TWK_BASE_MAP[256] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,3,0,0,0,2,0,0,0,0,0,0,4,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0};

const char TWK_BASE_MAP_INV[4] = {'A','T','G','C'};
const std::vector<std::string> TWK_SIMD_MAPPING = {"NONE","SSE","SSE2","SSE4","AVX","AVX2-256","AVX-512"};

inline void* aligned_malloc(size_t size, size_t align) {
    void *result;
    #ifdef _MSC_VER
    result = _aligned_malloc(size, align);
    #else
     if(posix_memalign(&result, align, size)) result = 0;
    #endif
    return result;
}

inline void aligned_free(void *ptr) {
    #ifdef _MSC_VER
        _aligned_free(ptr);
    #else
      free(ptr);
    #endif

}

// Unpacking types
#define TWK_LDD_NONE   0
#define TWK_LDD_VEC    1
#define TWK_LDD_LIST   2
#define TWK_LDD_BITMAP 4
#define TWK_LDD_ALL  ((TWK_LDD_VEC) | (TWK_LDD_LIST) | (TWK_LDD_BITMAP))

/****************************
*  SIMD definitions
****************************/
#if defined(_MSC_VER)
     /* Microsoft C/C++-compatible compiler */
     #include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
     /* GCC-compatible compiler, targeting x86/x86-64 */
     #include <x86intrin.h>
#elif defined(__GNUC__) && defined(__ARM_NEON__)
     /* GCC-compatible compiler, targeting ARM with NEON */
     #include <arm_neon.h>
#elif defined(__GNUC__) && defined(__IWMMXT__)
     /* GCC-compatible compiler, targeting ARM with WMMX */
     #include <mmintrin.h>
#elif (defined(__GNUC__) || defined(__xlC__)) && (defined(__VEC__) || defined(__ALTIVEC__))
     /* XLC or GCC-compatible compiler, targeting PowerPC with VMX/VSX */
     #include <altivec.h>
#elif defined(__GNUC__) && defined(__SPE__)
     /* GCC-compatible compiler, targeting PowerPC with SPE */
     #include <spe.h>
#endif

#if defined(__AVX512F__) && __AVX512F__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	6
#define SIMD_ALIGNMENT	64
#define GENOTYPE_TRIP_COUNT	256
#define SIMD_WIDTH      512
#elif defined(__AVX2__) && __AVX2__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	5
#define SIMD_ALIGNMENT	32
#define GENOTYPE_TRIP_COUNT	128
#define SIMD_WIDTH      256
#elif defined(__AVX__) && __AVX__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	4
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#define SIMD_WIDTH      128
#elif defined(__SSE4_1__) && __SSE4_1__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	3
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#define SIMD_WIDTH      128
#elif defined(__SSE2__) && __SSE2__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	2
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#define SIMD_WIDTH      128
#elif defined(__SSE__) && __SSE__ == 1
#define SIMD_AVAILABLE	0 // unsupported version
#define SIMD_VERSION	1
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#define SIMD_WIDTH      0
#else
#define SIMD_AVAILABLE	0
#define SIMD_VERSION	0
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#define SIMD_WIDTH      0
#endif

/****************************
*  Core genotype
****************************/
struct twk1_gt_t {
	twk1_gt_t() : n(0), miss(0){}
	virtual ~twk1_gt_t(){}
	virtual void print() =0;
	virtual void clear() =0;

	virtual uint32_t GetLength(const uint32_t p) const =0;
	virtual uint8_t GetRefByte(const uint32_t p) const =0;
	virtual uint8_t GetRefA(const uint32_t p) const =0;
	virtual uint8_t GetRefB(const uint32_t p) const =0;

	virtual twk1_gt_t* Clone() =0;
	virtual void Move(twk1_gt_t*& other) =0;
	virtual twk_buffer_t& AddBuffer(twk_buffer_t& buffer) const =0;
	virtual twk_buffer_t& ReadBuffer(twk_buffer_t& buffer) =0;

	friend inline twk_buffer_t& operator<<(twk_buffer_t& buffer, const twk1_gt_t& self){
		return(self.AddBuffer(buffer));
	}

	friend inline twk_buffer_t& operator>>(twk_buffer_t& buffer, twk1_gt_t& self){
		return(self.ReadBuffer(buffer));
	}

	uint32_t n: 31, miss: 1;
};

// internal genotype structure
template <class int_t>
struct twk1_igt_t : public twk1_gt_t {
	static const uint8_t psize = sizeof(int_t);

	twk1_igt_t() : data(nullptr){}
	~twk1_igt_t(){ delete[] data; }

	inline uint32_t GetLength(const uint32_t p) const{ return(data[p] >> (2 + 2*miss)); }
	inline uint8_t GetRefByte(const uint32_t p) const{ return(data[p] & ((1 << (2 + 2*miss)) - 1)); }
	inline uint8_t GetRefA(const uint32_t p) const{ return((data[p] >> (1 + miss)) & ((1 << (1 + miss)) - 1)); }
	inline uint8_t GetRefB(const uint32_t p) const{ return(data[p] & ((1 << (1 + miss)) - 1)); }

	twk_buffer_t& AddBuffer(twk_buffer_t& buffer) const{
		uint32_t n_write = this->n << 1 | this->miss;
		SerializePrimitive(n_write, buffer);
		for(int i = 0; i < this->n; ++i) SerializePrimitive<int_t>(this->data[i], buffer);
		return(buffer);
	}

	twk_buffer_t& ReadBuffer(twk_buffer_t& buffer){
		uint32_t n_write = 0;
		DeserializePrimitive(n_write, buffer);
		this->n = n_write >> 1;
		this->miss = n_write & 1;
		delete[] data; data = new int_t[n];
		for(int i = 0; i < this->n; ++i) DeserializePrimitive<int_t>(this->data[i], buffer);
		return(buffer);
	}

	/**<
	 * Clone operator for copying an inherited (virtual) class.
	 * @return
	 */
	twk1_gt_t* Clone(){
		twk1_igt_t* copied = new twk1_igt_t<int_t>;
		copied->n = n;
		copied->miss = miss;
		copied->data = new int_t[n];
		memcpy(copied->data, data, sizeof(int_t)*n);
		return(copied);
	}

	/**<
	 * Move operator for moving an inherited (virtual) class.
	 * @param other
	 */
	void Move(twk1_gt_t*& other){
		delete other;
		twk1_igt_t* dat = new twk1_igt_t<int_t>;
		//twk1_igt_t<int_t>* temp = reinterpret_cast<twk1_igt_t<int_t>*>(other);
		dat->n = n; n = 0;
		dat->miss = miss;
		std::swap(data, dat->data);
		other = dat;
		data = nullptr;
	}

	void clear(void){ delete[] data; data = nullptr; n = 0; }

	void print(){
		std::cerr << "print=" << n << std::endl;
		for(int i = 0; i < n; ++i){
			std::cerr << " <" << (data[i] & 3) << "," << (data[i] >> 2) << ">";
		}
		std::cerr << std::endl;
	}

	int_t* data;
};

/****************************
*  Core record
****************************/
struct twk1_t {
public:
	twk1_t();
	~twk1_t();

	twk1_t& operator=(const twk1_t& other);
	twk1_t& operator=(twk1_t&& other) noexcept;
	twk1_t(const twk1_t& other) = delete;
	twk1_t(twk1_t&& other) = delete;

	inline void EncodeAlleles(const char A, const char B){ this->alleles = (TWK_BASE_MAP[A] << 4) | TWK_BASE_MAP[B]; }
	inline char GetAlleleA() const{ return(TWK_BASE_MAP_INV[this->alleles >> 4]); }
	inline char GetAlleleB() const{ return(TWK_BASE_MAP_INV[this->alleles & 15]); }

	void clear();

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const twk1_t& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, twk1_t& self);

public:
    uint8_t  gt_ptype: 5, gt_flipped: 1, gt_phase: 1, gt_missing: 1;
    uint8_t  alleles;
    uint32_t pos, ac, an, rid, n_het, n_hom;
    twk1_gt_t* gt;
};

/****************************
*  Core containers
****************************/
struct twk1_block_t {
public:
	typedef twk1_block_t       self_type;
	typedef twk1_t             value_type;
	typedef value_type&        reference;
	typedef const value_type&  const_reference;
	typedef value_type*        pointer;
	typedef const value_type*  const_pointer;
	typedef std::ptrdiff_t     difference_type;
	typedef std::size_t        size_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
	twk1_block_t();
	twk1_block_t(const uint32_t p);
	~twk1_block_t();

	inline void operator+=(const twk1_t& rec){ this->Add(rec); }
	twk1_block_t& operator=(twk1_block_t&& other);
	twk1_block_t& operator=(const twk1_block_t& other) = delete;
	twk1_block_t(const twk1_block_t& other) = delete;
	twk1_block_t(twk1_block_t&& other) = delete;

	void Add(const twk1_t& rec);
	void resize(void);

	inline const uint32_t& size(void) const{ return(this->n); }

	inline reference front(void){ return(this->rcds[0]); }
	inline const_reference front(void) const{ return(this->rcds[0]); }
	inline reference back(void){ return(this->rcds[this->size() == 0 ? 0 : this->size() - 1]); }
	inline const_reference back(void) const{ return(this->rcds[this->size() == 0 ? 0 : this->size() - 1]); }

	inline reference operator[](const uint32_t position){ return(this->rcds[position]); }
	inline const_reference operator[](const uint32_t position) const{ return(this->rcds[position]); }
	inline reference at(const uint32_t position){ return(this->rcds[position]); }
	inline const_reference at(const uint32_t position) const{ return(this->rcds[position]); }

	inline pointer end(void){ return(&this->rcds[this->n]); }
	inline const_pointer end(void) const{ return(&this->rcds[this->n]); }

	void clear();

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const twk1_block_t& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, twk1_block_t& self);

public:
	uint32_t n, m, rid;
	uint32_t minpos, maxpos;
	twk1_t* rcds;
};

struct twk_oblock_t {
public:
	twk_oblock_t();
	void Write(std::ostream& stream, const uint32_t n, const uint32_t nc, const twk_buffer_t& buffer);

	friend std::ostream& operator<<(std::ostream& stream, const twk_oblock_t& self);
	friend std::istream& operator>>(std::istream& stream, twk_oblock_t& self);

public:
	uint32_t n, nc;
	twk_buffer_t bytes;
};

/****************************
*  Genotype specializations
****************************/
struct twk_igt_list {
public:
	struct ilist_t {
		ilist_t() :
			bin(0),
			vals(reinterpret_cast<uint64_t*>(aligned_malloc(2*sizeof(uint64_t), SIMD_ALIGNMENT)))
		{
			memset(vals, 0, sizeof(uint64_t)*2);
		}
		~ilist_t(){ aligned_free(vals); }

		inline void Set(uint32_t p){
			p %= 128;
			assert(p/64 < 2);
			vals[p/64] |= ((uint64_t)1 << (p % 64));
		}

		uint32_t bin;
		uint64_t* vals;
	};

	struct ilist_cont {
		ilist_cont(uint32_t n_s) : n(0), m(n_s*2), s(0), data(new ilist_t[n_s*2]){}
		ilist_cont() : n(0), m(0), s(0), data(nullptr){}
		~ilist_cont(){ delete[] data; }

		void resize(uint32_t n_s){
			n = 0; s = 0;
			delete[] data;
			m = std::ceil((float)2*n_s/128) + 1;
			data = new ilist_t[m];
			//std::cerr << "m=" << m << std::endl;
		}

		void resize2(uint32_t n_s){
			s = 0;
			delete[] data;
			m = n_s / 128 + 1;
			n = m;
			//std::cerr << "allocating=" << m << " registers" << std::endl;
			data = new ilist_t[m];
			//std::cerr << "m=" << m << std::endl;
		}

		void Add(const int32_t p){
			//std::cerr << "add=" << p << "/" << m*128 << " bin=" << p / 128 << "/" << m << std::endl;
			assert(p < m*128);
			//data[p / 128].Set(p);
			if(p / 128 == tail()){
				assert(n != 0);
				data[n - 1].Set(p);
			} else {
				//std::cerr << "p=" << p << std::endl;
				//std::cerr << "diff=" << p/128 << " and " << tail() << " n=" << n << "/" << m << std::endl;
				++n;
				assert(n < m);
				data[n - 1].bin = p / 128;
				data[n - 1].Set(p);
			}
			++s;
		}

		void operator+=(const int32_t p){
			//std::cerr << "add=" << p << "/" << m*128 << " bin=" << p / 128 << "/" << m << std::endl;
			assert(p < m*128);
			data[p / 128].Set(p);
			++s;
			if(pos.size()){
				if(pos.back() != p/128)
					pos.push_back(p/128);
			} else pos.push_back(p/128);
		}

		inline int32_t tail() const{ return(n == 0 ? -1 : data[n - 1].bin); }

		void cleanup(){
			if(s == 0){
				m = 0;
				n = 0;
				delete[] data;
				data = nullptr;
				pos.clear();
			}
		}

		uint32_t n, m, s;
		ilist_t* data;
		std::vector<uint32_t> pos;
	};

	struct ilist_cont_bins {
		ilist_cont_bins() : n(0), m(0), b(0), data(nullptr){}
		~ilist_cont_bins(){ delete[] data; }

		void resize(uint32_t n_s){
			delete[] data;
			m = 128;
			n = m;
			b = std::ceil((float)2*n_s / 128) + 1;
			//std::cerr << "registers=" << std::ceil((float)2*n_s / 128) + 1 << std::endl;
			b = std::ceil((float)b / m);
			data = new ilist_cont[m];
			//std::cerr << "allocating=" << 128 << " bins" << std::endl;
			for(int i = 0; i < m; ++i){
				//std::cerr << "resizing to=" << b*128 << std::endl;
				data[i].resize2(b*128);
			}
			//exit(1);
		}

		void operator+=(const int32_t p){
			//std::cerr << "Adding=" << p << "->" << p % (b * 128) << " in bucket " << p / (b * 128) << std::endl;
			data[p / (b * 128)] += p % (b * 128);
			if(pos.size()){
				if(pos.back() != (p / (b * 128)))
					pos.push_back(p / (b * 128));
			} else pos.push_back(p / (b * 128));
		}

		void cleanup(){
			for(int i = 0; i < m; ++i){
				data[i].cleanup();
			}
		}

		uint32_t n, m, b;
		//ilist_cont_bins* heap;
		ilist_cont* data; // containers for each bin
		std::vector<uint32_t> pos;
	};

	twk_igt_list() :
		own(1), n(0), m(0), l_list(0),
		list(nullptr), bv(nullptr), x_bins(0),
		bin_bitmap(reinterpret_cast<uint64_t*>(aligned_malloc(2*sizeof(uint64_t), SIMD_ALIGNMENT)))
	{
		memset(bin_bitmap, 0, sizeof(uint64_t)*2);
	}

	~twk_igt_list(){
		delete[] this->list;
		if(own) delete[] this->bv;
		aligned_free(bin_bitmap);
	}


	bool Build(const twk1_t& twk,
	           const uint32_t n_samples,
	           const bool resizeable = false)
		{
		if(n == 0){
			n = ceil((double)(n_samples*2)/64);
			if(resizeable){ m = twk.ac; }
			else { m = n_samples*2; }
			delete[] list; list = new uint32_t[m];
			if(own){
				bv = new uint64_t[n];
			}
			//std::cerr << "build now=" << n << "," << m << "," << l_list << std::endl;
		}
		//if(d.m == 0) d.resize(n_samples);

		if(twk.ac > m){
			//std::cerr << "resize=" << twk.ac << ">" << m << std::endl;
			m = twk.ac;
			delete[] list; list = new uint32_t[m];
		}

		if(own){
			memset(bv, 0, n*sizeof(uint64_t));
		}

		l_list = 0;
		uint32_t cumpos = 0;
		if(own){
			for(int i = 0; i < twk.gt->n; ++i){
				const uint32_t len  = twk.gt->GetLength(i);
				const uint8_t  refA = twk.gt->GetRefA(i);
				const uint8_t  refB = twk.gt->GetRefB(i);

				if(refA == 0 && refB == 0){
					cumpos += 2*len;
					continue;
				}

				for(int j = 0; j < 2*len; j+=2){
					if(refA != 0){ this->set(cumpos + j + 0); list[l_list++] = cumpos + j + 0; }
					if(refB != 0){ this->set(cumpos + j + 1); list[l_list++] = cumpos + j + 1; }
				}
				cumpos += 2*len;
			}
			assert(cumpos == n_samples*2);
			assert(l_list <= m);
		} else {
			//std::cerr << "in not own" << std::endl;

			//x.resize(n_samples);
			for(int i = 0; i < twk.gt->n; ++i){
				const uint32_t len  = twk.gt->GetLength(i);
				const uint8_t  refA = twk.gt->GetRefA(i);
				const uint8_t  refB = twk.gt->GetRefB(i);

				if(refA == 0 && refB == 0){
					cumpos += 2*len;
					continue;
				}

				for(int j = 0; j < 2*len; j+=2){
					if(refA == 1){
						list[l_list++] = cumpos + j + 0;
						//d.Add(cumpos + j + 0);
						//x += cumpos + j + 0;
					}
					if(refB == 1){
						list[l_list++] = cumpos + j + 1;
						//d.Add(cumpos + j + 1);
						//x += cumpos + j + 1;
					}
				}
				cumpos += 2*len;
			}
			assert(cumpos == n_samples*2);
			//std::cerr << "total bins=" << d.n << " with ac=" << twk.ac << std::endl;

			// Register positions.
			for(int i = 0; i < l_list; ++i){
				if(r_pos.size()){
					if(r_pos.back() != list[i] / 128)
						r_pos.push_back(list[i] / 128);
				} else r_pos.push_back(list[i] / 128);
			}

			// Convert register number into offset in 64-bit space.
			// This step is required to guarantee memory aligned lookups when
			// using SIMD instructions.
			for(int i = 0; i < r_pos.size(); ++i) r_pos[i] *= 128 / 64;

			//d.n = d.m;
			//std::cerr << "size=" << d.n << " and " << d.m << std::endl;
			//uint32_t divide = std::ceil((float)d.m / 2);
			uint32_t divide = std::ceil((float)2*n_samples/128) + 1;
			for(int i = 0; i < r_pos.size(); ++i){
				assert(r_pos[i] / divide < 2);
				bin_bitmap[r_pos[i]/divide] |= ((uint64_t)1 << ( (uint64_t)(((float)(r_pos[i] % divide) / divide)*64) ));
			}
			//std::cerr << std::bitset<64>(bin_bitmap[0]) << " " << std::bitset<64>(bin_bitmap[1]) << std::endl;

			/*for(int i = 0; i < d.n; ++i){
				aligned_free(d.data[i].vals);
				d.data[i].vals = nullptr;
			}*/


			/*x_bins = 0;
			for(int i = 0; i < x.n; ++i){
				x_bins += x.data[i].pos.size();
				//for(int j = 0; j < x.data[i].pos.size(); ++j){
					//std::cerr << i << "/" << x.n << " " << j << "/" << x.data[i].pos.size() << " -> " << x.data[i].pos[j] << std::endl;
				//}
			}
			//std::cerr << "done" << std::endl;
			//std::cerr << "bins=" << x_bins << " and " << x.pos.size() << std::endl;

			x.cleanup();*/
		}
		return true;
	}

	inline void reset(void){ memset(this->bv, 0, this->n*sizeof(uint64_t)); delete[] list; l_list = 0; }
	//inline const bool operator[](const uint32_t p) const{ return(this->bv[p] & (1L << (p % 8))); }
	/**<
	 * p >> 3 = p/8
	 * p % 8 = p & (8 - 1)
	 * @param p
	 * @return
	 */
	__attribute__((always_inline)) inline const bool get(const uint32_t p) const{ return(this->bv[(p >> 6)] & (1L << ( p & (64 - 1) )));}
	inline void set(const uint32_t p){ this->bv[p/64] |= (1L << (p % 64)); }
	//inline void set(const uint32_t p, const bool val){ this->bv[p/64] |= ((uint64_t)1 << (p % 64)); }

public:
	uint32_t own: 1, n: 31;
	uint32_t m, l_list; // n bytes, m allocated list entries
	uint32_t* list; // list entries
	uint64_t* bv; // bit-vector
	ilist_cont d;
	ilist_cont_bins x;
	uint32_t x_bins;
	// higher level bitmap checks
	uint64_t* bin_bitmap;
	std::vector<uint32_t> r_pos;
};

/**<
 * Data structure to representing a 1-bit allele representation of genotypes.
 * This data structure has to be aligned to `SIMD_ALIGNMENT` as specified
 * by the user CPU architecture. If no SIMD is available on the device then
 * use regular memory alignment.
 *
 * Special techniques to accelerate pairwise comparisons:
 * 1) Front and tail number of SIMD _elements_ (e.g. 128 bits / 16 bytes)
 *    that are either all 0 or 1. This allows the algorithm to either
 *    completely skip these stretches or resort to cheaper comparison functors.
 * 2) Counts of missingness needs to be maintained for these tail and head
 *    elements to function correctly.
 */
// INVERSE mask is cheaper in terms of instructions used.
// Lookup tables:
// 1 bit    2-bit
// 0|0 = 0  0|0 = 0
// 0|1 = 1  0|1 = 1
// 1|0 = 2  1|0 = 0100b = 4
// 1|1 = 3  1|1 = 0101b = 5
//          0|2 = 0010b = 2
//          1|2 = 0110b = 6
//          2|2 = 1010b = 12
//          2|0 = 1000b = 8
//          2|1 = 1001b = 9
const static uint8_t TWK_GT_BV_MASK_LOOKUP[16]      = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const static uint8_t TWK_GT_BV_MASK_MISS_LOOKUP[16] = {0, 0, 3, 3, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
const static uint8_t TWK_GT_BV_DATA_LOOKUP[16]      = {0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const static uint8_t TWK_GT_BV_DATA_MISS_LOOKUP[16] = {0, 1, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

struct twk_igt_vec {
public:
	twk_igt_vec():
		n(0),
		front_zero(0),
		tail_zero(0),
		data(nullptr),
		mask(nullptr)
	{
	}

	~twk_igt_vec(){
		aligned_free(data);
		aligned_free(mask);
	}

	inline void reset(void){ memset(this->data, 0, this->n*sizeof(uint64_t)); memset(this->mask, 0, this->n*sizeof(uint64_t)); }
	inline const bool operator[](const uint32_t p) const{ return(this->data[p] & (1L << (p % 64))); }
	inline const bool get(const uint32_t p) const{ return(this->data[p/64] & (1L << (p % 64)));}
	inline void SetData(const uint32_t p){ this->data[p/64] |= (1L << (p % 64)); }
	inline void SetData(const uint32_t p, const bool val){ this->data[p/64] |= ((uint64_t)val << (p % 64)); }
	inline void SetMask(const uint32_t p){ this->mask[p/64] |= (1L << (p % 64)); }
	inline void SetMask(const uint32_t p, const bool val){ this->mask[p/64] |= ((uint64_t)val << (p % 64)); }

	/**<
	 *
	 * @param rec
	 * @param n_samples
	 * @param alignment
	 * @return
	 */
	bool Build(const twk1_t& rec,
	           const uint32_t& n_samples)
	{
		if(n == 0){
			n = ceil((double)(n_samples*2)/64);
			n += (n*64) % 128; // must be divisible by 128-bit register
			data = reinterpret_cast<uint64_t*>(aligned_malloc(n*sizeof(uint64_t), SIMD_ALIGNMENT));
			if(rec.gt_missing){
				mask = reinterpret_cast<uint64_t*>(aligned_malloc(n*sizeof(uint64_t), SIMD_ALIGNMENT));
			}
		}

		memset(data, 0, n*sizeof(uint64_t));
		if(rec.gt_missing)
			memset(mask, 0, n*sizeof(uint64_t));

		uint32_t cumpos = 0;
		for(int i = 0; i < rec.gt->n; ++i){
			const uint32_t len  = rec.gt->GetLength(i);
			const uint8_t  refA = rec.gt->GetRefA(i);
			const uint8_t  refB = rec.gt->GetRefB(i);

			if(refA == 0 && refB == 0){
				cumpos += 2*len;
				continue;
			}

			for(int j = 0; j < 2*len; j+=2){
				if(refA == 1){ this->SetData(cumpos + j + 0); }
				if(refB == 1){ this->SetData(cumpos + j + 1); }
				if(refA == 2){ this->SetMask(cumpos + j + 0); this->SetMask(cumpos + j + 1); }
				if(refB == 2){ this->SetMask(cumpos + j + 0); this->SetMask(cumpos + j + 1); }
			}
			cumpos += 2*len;
		}

		if(rec.gt_missing){
			// Invert mask
			/*for(int i = 0; i < n; ++i){
				mask[i] = ~mask[i];
			}*/
		}
		assert(cumpos == n_samples*2);


		const uint32_t byteAlignedEnd  = (2*n_samples/SIMD_WIDTH) * (SIMD_WIDTH / 64);

		int j = 0;
		if(rec.gt_missing){
			// Search from left->right
			for(; j < byteAlignedEnd; ++j){
				if(this->data[j] != 0 || this->mask[j] != 0)
					break;
			}

			// Front of zeroes
			this->front_zero = ((j - 1 < 0 ? 0 : j - 1)*4)/GENOTYPE_TRIP_COUNT;
			if(j != byteAlignedEnd){
				j = byteAlignedEnd - 1;
				for(; j > 0; --j){
					if(this->data[j] != 0 || this->mask[j] != 0)
						break;
				}
			}
		}
		// If no data is missing then we have no mask. Run this subroutine instead
		// to avoid memory problems.
		else {
			// Search from left->right
			for(; j < byteAlignedEnd; ++j){
				if(this->data[j] != 0)
					break;
			}

			// Front of zeroes
			this->front_zero = ((j - 1 < 0 ? 0 : j - 1)*4)/GENOTYPE_TRIP_COUNT;
			if(j != byteAlignedEnd){
				j = byteAlignedEnd - 1;
				for(; j > 0; --j){
					if(this->data[j] != 0)
						break;
				}
			}
		}

		// Tail of zeroes
		this->tail_zero = byteAlignedEnd == j ? 0 : ((byteAlignedEnd - (j+1))*4)/GENOTYPE_TRIP_COUNT;

		return(true);
	}

public:
	uint32_t  n; // n bytes, m allocated entries
	uint16_t front_zero; // leading zeros in aligned vector width
	uint16_t tail_zero;  // trailing zeros in aligned vector width
	uint64_t* data;
	uint64_t* mask;
};

// gtocc
struct twk_gtocc {
	twk_gtocc(): n(0), cumpos(nullptr){}
	twk_gtocc(const uint32_t n_s): n(n_s+1), cumpos(new uint32_t*[n]){
		for(int i = 0; i < n; ++i){
			cumpos[i] = new uint32_t[16];
			memset(cumpos[i], 0, sizeof(uint32_t)*16);
		}
	}
	~twk_gtocc(){
		for(int i = 0; i < n; ++i) delete [] cumpos[i];
		delete[] cumpos;
	}

	void operator=(const twk1_gt_t* gt){
		assert(n != 0);
		for(int i = 0; i < n; ++i) memset(cumpos[i], 0, sizeof(uint32_t)*16);
		uint32_t cum = 0;
		for(int i = 0; i < gt->n; ++i){
			const uint32_t l = gt->GetLength(i);
			const uint8_t ref = gt->GetRefByte(i);
			for(int j = 0; j < l; ++j){

			}
		}
	}

	uint32_t n;
	// <0,0> <0,1> <0,2> <1,0> <1,1> <1,2> <2,0> <2,1> <2,2>
	uint32_t** cumpos;
};


// output
struct twk1_two_t {
public:
	const static uint32_t packed_size = sizeof(uint16_t) + 2*sizeof(int32_t) + 2*sizeof(uint32_t) + 11*sizeof(double);

	twk1_two_t() :
		controller(0), ridA(0), ridB(0), Amiss(0), Aphased(0), Apos(0), Bmiss(0), Bphased(0), Bpos(0), R(0), R2(0), D(0), Dprime(0), P(0), ChiSqModel(0), ChiSqFisher(0){ memset(cnt, 0, sizeof(double)*4); }
	~twk1_two_t() = default;

	inline double& operator[](const uint32_t& p){ return(this->cnt[p]); }
	inline const double& operator[](const uint32_t& p) const{ return(this->cnt[p]); }

	void PrintUnphasedCounts(void) const;
	void PrintPhasedCounts(void) const;

	friend twk_buffer_t& operator<<(twk_buffer_t& os, const twk1_two_t& entry){
		SerializePrimitive(entry.controller, os);
		SerializePrimitive(entry.ridA, os);
		SerializePrimitive(entry.ridB, os);
		uint32_t packA = entry.Apos << 2 | entry.Aphased << 1 | entry.Amiss;
		uint32_t packB = entry.Bpos << 2 | entry.Bphased << 1 | entry.Bmiss;
		SerializePrimitive(packA, os);
		SerializePrimitive(packB, os);
		SerializePrimitive(entry.cnt[0], os);
		SerializePrimitive(entry.cnt[1], os);
		SerializePrimitive(entry.cnt[2], os);
		SerializePrimitive(entry.cnt[3], os);
		SerializePrimitive(entry.D, os);
		SerializePrimitive(entry.Dprime, os);
		SerializePrimitive(entry.R, os);
		SerializePrimitive(entry.R2, os);
		SerializePrimitive(entry.P, os);
		SerializePrimitive(entry.ChiSqFisher, os);
		SerializePrimitive(entry.ChiSqModel, os);
		return os;
	}

	friend twk_buffer_t& operator>>(twk_buffer_t& os, twk1_two_t& entry){
		DeserializePrimitive(entry.controller, os);
		DeserializePrimitive(entry.ridA, os);
		DeserializePrimitive(entry.ridB, os);
		uint32_t packA = 0;
		uint32_t packB = 0;
		DeserializePrimitive(packA, os);
		DeserializePrimitive(packB, os);
		entry.Apos = packA >> 2;
		entry.Bpos = packB >> 2;
		entry.Aphased = (packA >> 1) & 1;
		entry.Bphased = (packB >> 1) & 1;
		entry.Amiss = packA & 1;
		entry.Bmiss = packB & 1;
		DeserializePrimitive(entry.cnt[0], os);
		DeserializePrimitive(entry.cnt[1], os);
		DeserializePrimitive(entry.cnt[2], os);
		DeserializePrimitive(entry.cnt[3], os);
		DeserializePrimitive(entry.D, os);
		DeserializePrimitive(entry.Dprime, os);
		DeserializePrimitive(entry.R, os);
		DeserializePrimitive(entry.R2, os);
		DeserializePrimitive(entry.P, os);
		DeserializePrimitive(entry.ChiSqFisher, os);
		DeserializePrimitive(entry.ChiSqModel, os);
		return(os);
	}

	/**<
	 * Print out record in uncompressed LD format.
	 * @param os
	 * @return
	 */
	std::ostream& PrintLD(std::ostream& os) const {
		os << controller << '\t' << ridA << '\t' << Apos << '\t' << ridB << '\t' << Bpos << '\t'
		   << cnt[0] << '\t' << cnt[1] << '\t' << cnt[2] << '\t' << cnt[3] << '\t' << D << '\t'
		   << Dprime << '\t' << R << '\t' << R2 << '\t' << P << '\t' << ChiSqFisher << '\t' << ChiSqModel << '\n';
		return(os);
	}

	/**<
	 * Print out JSON
	 * @param os
	 * @return
	 */
	std::ostream& PrintLDJson(std::ostream& os) const {
		os << '[' << controller << ',' << ridA << ',' << Apos << ',' << ridB << ',' << Bpos << ','
		   << cnt[0] << ',' << cnt[1] << ',' << cnt[2] << ',' << cnt[3] << ',' << D << ','
		   << Dprime << ',' << R << ',' << R2 << ',' << P << ',' << ChiSqFisher << ',' << ChiSqModel << ']' << '\n';
		return(os);
	}

	friend std::ostream& operator<<(std::ostream& os, const twk1_two_t& entry){
		os << entry.controller << '\t' << entry.ridA << '\t' << entry.Apos << '\t' << entry.ridB << '\t' << entry.Bpos << '\t'
		   << entry.cnt[0] << '\t' << entry.cnt[1] << '\t' << entry.cnt[2] << '\t' << entry.cnt[3] << '\t'
		   << entry.D << '\t' << entry.Dprime << '\t' << entry.R << '\t' << entry.R2 << '\t' << entry.P << '\t' << entry.ChiSqFisher << '\t' << entry.ChiSqModel;
		return os;
	}

	inline void SetUsedPhasedMath(const bool yes = true)   { this->controller |= yes << 0;  } // 1
	inline void SetSameContig(const bool yes = true)       { this->controller |= yes << 1;  } // 2
	inline void SetLongRange(const bool yes = true)        { this->controller |= yes << 2;  } // 4
	inline void SetCompleteLD(const bool yes = true)       { this->controller |= yes << 3;  } // 8
	inline void SetPerfectLD(const bool yes = true)        { this->controller |= yes << 4;  } // 16
	inline void SetMultipleRoots(const bool yes = true)    { this->controller |= yes << 5;  } // 32
	inline void SetFastMode(const bool yes = true)         { this->controller |= yes << 6;  } // 64
	inline void SetSampled(const bool yes = true)          { this->controller |= yes << 7;  } // 128
	inline void SetHasMissingValuesA(const bool yes = true){ this->controller |= yes << 8;  } // 256
	inline void SetHasMissingValuesB(const bool yes = true){ this->controller |= yes << 9;  } // 512
	inline void SetLowACA(const bool yes = true)           { this->controller |= yes << 10; } // 1024
	inline void SetLowACB(const bool yes = true)           { this->controller |= yes << 11; } // 2048

	void clear(){
		controller = 0; ridA = 0; ridB = 0;
		Amiss = 0; Aphased = 0; Apos = 0;
		Bmiss = 0; Bphased = 0; Bpos = 0;
		R = 0; R2 = 0; D = 0; Dprime = 0; P = 0;
		ChiSqModel = 0; ChiSqFisher = 0;
		memset(cnt, 0, sizeof(double)*4);
	}

	bool operator<(const twk1_two_t& other) const{
		if(ridA < other.ridA) return true;
		if(other.ridA < ridA) return false;
		if(ridB < other.ridB) return true;
		if(other.ridB < ridB) return false;
		if(Apos < other.Apos) return true;
		if(other.Apos < Apos) return false;
		if(Bpos < other.Bpos) return true;
		if(other.Bpos < Bpos) return false;
		return false;
	}

public:
	uint16_t controller;
	uint32_t ridA, ridB;
	uint32_t Amiss: 1, Aphased: 1, Apos: 30;
	uint32_t Bmiss: 1, Bphased: 1, Bpos: 30;
	double R, R2, D, Dprime, P;
	double ChiSqModel;   // Chi-Squared critical value for 3x3 contingency table
	double ChiSqFisher;  // Chi-Squared critical value for 2x2 contingency table
	double cnt[4];
};


struct twk_oblock_two_t {
public:
	twk_oblock_two_t() : n(0), nc(0){}

	inline void operator+=(const twk1_two_t& entry){ bytes << entry; }

	void Write(std::ostream& stream, const uint32_t n, const uint32_t nc, const twk_buffer_t& buffer){
		uint8_t marker = 1;
		SerializePrimitive(marker, stream);
		SerializePrimitive(n, stream);
		SerializePrimitive(nc, stream);
		stream.write(buffer.data(), buffer.size());
	}

	friend std::ostream& operator<<(std::ostream& stream, const twk_oblock_two_t& self){
		uint8_t marker = 1;
		SerializePrimitive(marker, stream);
		SerializePrimitive(self.n, stream);
		SerializePrimitive(self.nc, stream);
		assert(self.nc == self.bytes.size());
		stream.write(self.bytes.data(), self.bytes.size());
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, twk_oblock_two_t& self){
		// marker has to be read outside
		self.bytes.reset();
		DeserializePrimitive(self.n,  stream);
		DeserializePrimitive(self.nc, stream);
		self.bytes.reset(); self.bytes.resize(self.nc);
		stream.read(self.bytes.data(), self.nc);
		self.bytes.n_chars_ = self.nc;
		return(stream);
	}

public:
	uint32_t n, nc;
	twk_buffer_t bytes;
};

struct twk1_two_block_t {
public:
	typedef twk1_two_block_t   self_type;
	typedef twk1_two_t         value_type;
	typedef value_type&        reference;
	typedef const value_type&  const_reference;
	typedef value_type*        pointer;
	typedef const value_type*  const_pointer;
	typedef std::ptrdiff_t     difference_type;
	typedef std::size_t        size_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
   	twk1_two_block_t() : n(0), m(0), rcds(nullptr){}
   	twk1_two_block_t(const uint32_t p): n(0), m(p), rcds(new twk1_two_t[p]){}
	~twk1_two_block_t(){ delete[] rcds; }

	inline twk1_two_block_t& operator+=(const twk1_two_t& rec){ return(this->Add(rec)); }
	twk1_two_block_t& Add(const twk1_two_t& rec){
		if(n == m) this->resize(); // will also trigger when n=0 and m=0
		rcds[n++] = rec;
		return(*this);
	}

	void resize(const uint32_t p){
		if(p < n) return;

		//std::cerr << "resizing to " << p << " @ " << n << "/" << m << std::endl;

		twk1_two_t* temp = rcds;
		rcds = new twk1_two_t[p];
		for(int i = 0; i < n; ++i) rcds[i] = std::move(temp[i]); // move records over
		delete[] temp;
		m = p;
	}

	void resize(void){
		if(this->rcds == nullptr){
			this->rcds = new twk1_two_t[500];
			this->n = 0;
			this->m = 500;
			return;
		}
		//std::cerr << "resizing=" << n << "/" << m << std::endl;

		twk1_two_t* temp = rcds;
		rcds = new twk1_two_t[m*2];
		for(int i = 0; i < n; ++i) rcds[i] = std::move(temp[i]); // move records over
		delete[] temp;
		m*=2;
	}

	void reserve(const uint32_t p);

	inline const uint32_t& size(void) const{ return(this->n); }

	inline reference front(void){ return(this->rcds[0]); }
	inline const_reference front(void) const{ return(this->rcds[0]); }
	inline reference back(void){ return(this->rcds[this->size() == 0 ? 0 : this->size() - 1]); }
	inline const_reference back(void) const{ return(this->rcds[this->size() == 0 ? 0 : this->size() - 1]); }

	inline reference operator[](const uint32_t position){ return(this->rcds[position]); }
	inline const_reference operator[](const uint32_t position) const{ return(this->rcds[position]); }
	inline reference at(const uint32_t position){ return(this->rcds[position]); }
	inline const_reference at(const uint32_t position) const{ return(this->rcds[position]); }

	inline pointer start(void){ return(&this->rcds[0]); }
	inline const_pointer start(void) const{ return(&this->rcds[0]); }
	inline pointer end(void){ return(&this->rcds[this->n]); }
	inline const_pointer end(void) const{ return(&this->rcds[this->n]); }


	void reset(){
		n = 0;
	}

	void clear(){
		delete[] rcds; rcds = nullptr;
		n = 0; m = 0;
	}

	bool Sort(){
		//std::cerr << "sorting=" << n << std::endl;
		std::sort(start(), end());
		return(true);
	}

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const twk1_two_block_t& self){
		SerializePrimitive(self.n, buffer);
		SerializePrimitive(self.m, buffer);
		for(int i = 0; i < self.n; ++i) buffer << self.rcds[i];
		return(buffer);
	}

	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, twk1_two_block_t& self){
		delete[] self.rcds; self.rcds = nullptr;
		DeserializePrimitive(self.n, buffer);
		DeserializePrimitive(self.m, buffer);
		self.rcds = new twk1_two_t[self.m];
		for(int i = 0; i < self.n; ++i) buffer >> self.rcds[i];
		return(buffer);
	}

public:
	uint32_t n, m;
	twk1_two_t* rcds;
};

/**<
 * Settings/parameters for both `twk_ld` and `twk_ld_engine`. Placed in this header
 * to circumvent linkage problems in headers.
 */
struct twk_ld_settings {
public:
	twk_ld_settings();
	std::string GetString() const;

public:
	bool square, window, low_memory, bitmaps; // using square compute, using window compute
	bool force_phased, forced_unphased, force_cross_intervals;
	int32_t c_level, bl_size, b_size, l_window; // compression level, block_size, output block size, window size in bp
	int32_t n_threads, cycle_threshold, ldd_load_type;
	std::string in, out; // input file, output file/cout
	double minP, minR2, maxR2, minDprime, maxDprime;
	int32_t n_chunks, c_chunk;
	std::vector<std::string> ival_strings; // unparsed interval strings
};

}

#endif /* TWK_CORE_H_ */
