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
#include <atomic>

#include "tomahawk.h"
#include "buffer.h"
#include "header.h"

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

class SpinLock{
public:
    void lock(){
        while(lck.test_and_set(std::memory_order_acquire))
        {}
    }
    void unlock(){lck.clear(std::memory_order_release);}

private:
    std::atomic_flag lck = ATOMIC_FLAG_INIT;
};

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

	/**<
	 * 	This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
	 * 	Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
	 * 	Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
	 *
	 * 	Originally code written by Jan Wigginton
	 * 	Modified to use Tomahawk data by Marcus D. R. Klarqvist
	 */
	void calculateHardyWeinberg(void);

	void clear();

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const twk1_t& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, twk1_t& self);

public:
    uint8_t  gt_ptype: 5, gt_flipped: 1, gt_phase: 1, gt_missing: 1;
    uint8_t  alleles;
    uint32_t pos, ac, an, rid, n_het, n_hom;
    double hwe;
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
		own(1), n(0), m(0), l_list(0), l_aa(0), l_het(0),
		list(nullptr), list_aa(nullptr), list_het(nullptr), bv(nullptr),
		x_bins(0),
		bin_bitmap(reinterpret_cast<uint64_t*>(aligned_malloc(2*sizeof(uint64_t), SIMD_ALIGNMENT)))
	{
		memset(bin_bitmap, 0, sizeof(uint64_t)*2);
	}

	~twk_igt_list(){
		delete[] this->list;
		delete[] this->list_aa;
		delete[] this->list_het;
		if(own) delete[] this->bv;
		aligned_free(bin_bitmap);
	}


	bool Build(const twk1_t& twk,
	           const uint32_t n_samples,
	           const bool resizeable = false,
			   bool build_unphased = true)
		{
		if(n == 0){
			n = ceil((double)(n_samples*2)/64);
			if(resizeable){ m = twk.ac; }
			else { m = n_samples*2; }
			delete[] list; list = new uint32_t[m];
			if(build_unphased){
				delete[] list_aa; list_aa = new uint32_t[m];
				delete[] list_het; list_het = new uint32_t[m];
			}
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
			delete[] list_aa; list_aa = nullptr;
			delete[] list_het; list_het = nullptr;

			if(build_unphased){
				list_aa = new uint32_t[m];
				list_het = new uint32_t[m];
			}
		}

		if(own){
			memset(bv, 0, n*sizeof(uint64_t));
		}

		l_list = 0; l_aa = 0; l_het = 0;
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
					if(build_unphased && refA == 1 && refB == 1){ list_aa[l_aa++] = cumpos + j + 0; }
					if(build_unphased && ((refA == 0 && refB == 1) || (refA == 1 && refB == 0))){ list_het[l_het++] = cumpos + j + 0; }
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
					if(refA != 0){ this->set(cumpos + j + 0); list[l_list++] = cumpos + j + 0; }
					if(refB != 0){ this->set(cumpos + j + 1); list[l_list++] = cumpos + j + 1; }
					if(build_unphased && refA == 1 && refB == 1){ list_aa[l_aa++] = cumpos + j + 0; }
					if(build_unphased && ((refA == 0 && refB == 1) || (refA == 1 && refB == 0))){ list_het[l_het++] = cumpos + j + 0; }
				}
				cumpos += 2*len;
			}
			assert(cumpos == n_samples*2);
			//std::cerr << "total bins=" << d.n << " with ac=" << twk.ac << std::endl;

			// Register positions.
			r_pos.clear();
			for(int i = 0; i < l_list; ++i){
				if(r_pos.size()){
					if(r_pos.back() != list[i] / 128)
						r_pos.push_back(list[i] / 128);
				} else r_pos.push_back(list[i] / 128);
			}

			r_aa.clear(); r_het.clear();
			if(build_unphased){
				for(int i = 0; i < l_aa; ++i){
					if(r_aa.size()){
						if(r_aa.back() != list_aa[i] / 128)
							r_aa.push_back(list_aa[i] / 128);
					} else r_aa.push_back(list_aa[i] / 128);
				}

				for(int i = 0; i < l_het; ++i){
					if(r_het.size()){
						if(r_het.back() != list_het[i] / 128)
							r_het.push_back(list_het[i] / 128);
					} else r_het.push_back(list_het[i] / 128);
				}
			}
			//std::cerr << l_list << "," << l_het << "," << l_aa << " and " << r_pos.size() << "," << r_aa.size() << "," << r_het.size() << std::endl;


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

		// Delete list as we use registers
		delete[] list; list = nullptr;
		delete[] list_aa; list_aa = nullptr;
		delete[] list_het; list_het = nullptr;
		//l_list = 0;


		return true;
	}

	inline void reset(void){ memset(this->bv, 0, this->n*sizeof(uint64_t)); delete[] list; l_list = 0; }
	//inline const bool operator[](const uint32_t p) const{ return(this->bv[p] & (1L << (p % 8))); }

	__attribute__((always_inline)) inline const bool get(const uint32_t p) const{ return(this->bv[(p >> 6)] & (1L << ( p & (64 - 1) )));}
	inline void set(const uint32_t p){ this->bv[p/64] |= (1L << (p % 64)); }
	//inline void set(const uint32_t p, const bool val){ this->bv[p/64] |= ((uint64_t)1 << (p % 64)); }

public:
	uint32_t own: 1, n: 31;
	uint32_t m, l_list, l_aa, l_het; // n bytes, m allocated list entries
	uint32_t* list, *list_aa, *list_het; // list entries
	uint64_t* bv; // bit-vector
	//ilist_cont d;
	//ilist_cont_bins x;
	uint32_t x_bins;
	// higher level bitmap checks
	uint64_t* bin_bitmap;
	std::vector<uint32_t> r_pos, r_aa, r_het;
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
	twk_igt_vec();
	~twk_igt_vec();

	inline void reset(void){ memset(this->data, 0, this->n*sizeof(uint64_t)); memset(this->mask, 0, this->n*sizeof(uint64_t)); }
	inline const bool operator[](const uint32_t p) const{ return(this->data[p] & (1L << (p % 64))); }
	inline const bool get(const uint32_t p) const{ return(this->data[p/64] & (1L << (p % 64)));}
	inline void SetData(const uint32_t p){ this->data[p/64] |= (1L << (p % 64)); }
	inline void SetData(const uint32_t p, const bool val){ this->data[p/64] |= ((uint64_t)val << (p % 64)); }
	inline void SetMask(const uint32_t p){ this->mask[p/64] |= (1L << (p % 64)); }
	inline void SetMask(const uint32_t p, const bool val){ this->mask[p/64] |= ((uint64_t)val << (p % 64)); }

	/**<
	 * Constructs a 1-bitvector from from the given compressed genotypes. The
	 * bitvector is guaranteed to fit in the largest SIMD vector available and
	 * has enforced memory boundary alignments according to the hardware.
	 * @param rec       Input reference twk1_t record
	 * @param n_samples Total number of samples.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	bool Build(const twk1_t& rec, const uint32_t& n_samples);

public:
	uint32_t  n; // n bytes, m allocated entries
	uint16_t front_zero; // leading zeros in aligned vector width
	uint16_t tail_zero;  // trailing zeros in aligned vector width
	uint64_t* data;
	uint64_t* mask;
};

// output
struct twk1_two_t {
public:
	const static uint32_t packed_size = sizeof(uint16_t) + 2*sizeof(int32_t) +
	                                    2*sizeof(uint32_t) + 11*sizeof(double);

	twk1_two_t();
	~twk1_two_t() = default;

	inline double& operator[](const uint32_t& p){ return(this->cnt[p]); }
	inline const double& operator[](const uint32_t& p) const{ return(this->cnt[p]); }

	void PrintUnphasedCounts(void) const;
	void PrintPhasedCounts(void) const;

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
	inline void SetInvalidHWEA(const bool yes = true)      { this->controller |= yes << 12; } // 4096
	inline void SetInvalidHWEB(const bool yes = true)      { this->controller |= yes << 13; } // 8192

	void clear();
	bool operator<(const twk1_two_t& other) const;

	/**<
	 * Print out record in uncompressed LD format. Note that contig
	 * identifiers will be written in lieu of contig names.
	 * @param os Ouput ostream reference.
	 * @return   Returns the output ostream reference.
	 */
	std::ostream& PrintLD(std::ostream& os, VcfHeader* hdr) const;

	/**<
	 * Print out JSON to a given ostream. Note that contig identifiers
	 * will be written in lieu of contig names.
	 * @param os Ouput ostream reference.
	 * @return   Returns the output ostream reference.
	 */
	std::ostream& PrintLDJson(std::ostream& os) const;

	friend twk_buffer_t& operator<<(twk_buffer_t& os, const twk1_two_t& entry);
	friend twk_buffer_t& operator>>(twk_buffer_t& os, twk1_two_t& entry);
	friend std::ostream& operator<<(std::ostream& os, const twk1_two_t& entry);

public:
	/**<
	 * The controller is a 16-bit vector of flags:
	 *    1: Used phased math.
	 *    2: Acceptor and donor variants are on the same contig.
	 *    3: Acceptor and donor variants are far apart on the same contig.
	 *    4: The output contingency matrix has at least one empty cell (referred to as complete).
	 *    5: Output correlation coefficient is perfect (1.0).
	 *    6: Output solution is one of >1 possible solutions. This only occurs for unphased pairs.
	 *    7: Output data was generated in 'fast mode'.
	 *    8: Output data is estimated from a subsampling of the total pool of genotypes.
	 *    9: Donor vector has missing value(s).
	 *    10: Acceptor vector has missing value(s).
	 *    11: Donor vector has low allele count (<5).
	 *    12: Acceptor vector has low allele count (<5).
	 *    13: Acceptor vector has a HWE-P value < 1e-4.
	 *    14: Donor vector has a HWE-P value < 1e-4.
	 */
	uint16_t controller;
	uint32_t ridA, ridB;
	uint32_t Amiss: 1, Aphased: 1, Apos: 30;
	uint32_t Bmiss: 1, Bphased: 1, Bpos: 30;
	double R, R2, D, Dprime, P;
	double ChiSqModel;   // Chi-Squared critical value for 3x3 contingency table
	double ChiSqFisher;  // Chi-Squared critical value for 2x2 contingency table
	double cnt[4];       // REFREF,REFALT,ALTREF,ALTALT
};


struct twk_oblock_two_t {
public:
	twk_oblock_two_t();

	inline void operator+=(const twk1_two_t& entry){ bytes << entry; }

	void Write(std::ostream& stream, const uint32_t n, const uint32_t nc, const twk_buffer_t& buffer);

	friend std::ostream& operator<<(std::ostream& stream, const twk_oblock_two_t& self);
	friend std::istream& operator>>(std::istream& stream, twk_oblock_two_t& self);

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

public:
   	twk1_two_block_t();
   	twk1_two_block_t(const uint32_t p);
	~twk1_two_block_t();

	inline twk1_two_block_t& operator+=(const twk1_two_t& rec){ return(this->Add(rec)); }
	twk1_two_block_t& Add(const twk1_two_t& rec);

	void resize(const uint32_t p);
	void resize(void);
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

	void reset();
	void clear();
	bool Sort();

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const twk1_two_block_t& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, twk1_two_block_t& self);

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
	bool square, window, low_memory, bitmaps, single; // using square compute, using window compute
	bool force_phased, forced_unphased, force_cross_intervals;
	int32_t c_level, bl_size, b_size, l_window; // compression level, block_size, output block size, window size in bp
	int32_t n_threads, cycle_threshold, ldd_load_type;
	int32_t l_surrounding; // left,right-padding in base-pairs when running in single mode
	std::string in, out; // input file, output file/cout
	double minP, minR2, maxR2, minDprime, maxDprime;
	int32_t n_chunks, c_chunk;
	std::vector<std::string> ival_strings; // unparsed interval strings
};

/**<
 * Standard summary statistics object.
 */
struct twk_sstats {
    // Functional pointer definitions used in reduce/aggregate subroutines.
    typedef void (twk_sstats::*aggfunc)(const twk1_two_t*);
    typedef double (twk_sstats::*redfunc)(const uint32_t) const;

    twk_sstats() : n(0), total(0), total_squared(0), min(0), max(0){}

    template <class T> void Add(const T value, const double weight = 1){
        this->total         += value;
        this->total_squared += value*value;
        this->n             += weight;
        this->min = value < min ? value : min;
        this->max = value > max ? value : max;
    }

    void AddR2(const twk1_two_t* rec)  { Add(rec->R2); }
    void AddR(const twk1_two_t* rec)   { Add(rec->R);  }
    void AddD(const twk1_two_t* rec)   { Add(rec->D);  }
    void AddDprime(const twk1_two_t* rec){ Add(rec->Dprime); }
    void AddP(const twk1_two_t* rec)   { Add(rec->P);  }
    void AddHets(const twk1_two_t* rec){
        Add((rec->cnt[1] + rec->cnt[2]) / (rec->cnt[0] + rec->cnt[1] + rec->cnt[2] + rec->cnt[3]));
    }

    void AddAlts(const twk1_two_t* rec){
        Add((rec->cnt[3]) / (rec->cnt[0] + rec->cnt[1] + rec->cnt[2] + rec->cnt[3]));
    }

    double GetMean(const uint32_t min = 0) const {
        if(n < min || min == 0) return(0);
        return(total / n);
    }

    double GetCount(const uint32_t min = 0) const {
        if(n < min) return(0);
        return(n);
    }

    double GetStandardDeviation(const uint32_t cutoff = 2) const{
        if(this->n < cutoff) return(0);
        return(sqrt(this->total_squared/this->n - (this->total / this->n)*(this->total / this->n)));
    }

    // Accessor functions
    inline double GetTotal(const uint32_t cutoff = 0) const{ return(total < cutoff ? 0 : total); }
    inline double GetTotalSquared(const uint32_t cutoff = 0) const{ return(total_squared < cutoff ? 0 : total_squared); }
    inline double GetMin(const uint32_t cutoff = 0) const{ return(this->min); }
    inline double GetMax(const uint32_t cutoff = 0) const{ return(this->max); }

    void operator+=(const twk_sstats& other){
        n += other.n;
        total += other.total;
        total_squared += other.total_squared;
        min = std::min(min, other.min);
        max = std::max(max, other.max);
    }

public:
    uint64_t n;
    double total, total_squared;
    double min, max;
};

/**<
 * Structure for operating on aggregated datasets.
 */
struct twk1_aggregate_t {
public:
	struct offset_tuple {
	    offset_tuple() : range(0), min(std::numeric_limits<uint32_t>::max()), max(0){}
	    uint64_t range;
	    uint32_t min, max;
	};

public:
	twk1_aggregate_t();
	twk1_aggregate_t(const uint32_t x, const uint32_t y);
	~twk1_aggregate_t();

	friend std::ostream& operator<<(std::ostream& stream, const twk1_aggregate_t& agg);

	bool Open(std::string input);

public:
	// magic header
	uint32_t n, x, y, bpx, bpy, n_original;
	uint64_t range;
	std::string filename; // input filename
	std::vector<offset_tuple> rid_offsets; // mat offsets
	double* data;
	// EOF
};

}

#endif /* TWK_CORE_H_ */
