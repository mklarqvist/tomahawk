#ifndef TWK_LD_H_
#define TWK_LD_H_

#include <cassert>
#include <thread>

#include "core.h"
#include "twk_reader.h"
#include "fisher_math.h"
#include "spinlock.h"
#include "timer.h"
#include "writer.h"
#include "intervals.h"

// Method 0: None
// Method 1: Performance + no math
// Method 2: Debug correctness
// Method 3: Print LD difference Phased and Unphased
#define SLAVE_DEBUG_MODE 0

namespace tomahawk {

#define TWK_LD_SIMD_REFREF 0
#define TWK_LD_SIMD_ALTREF 1
#define TWK_LD_SIMD_REFALT 2
#define TWK_LD_SIMD_ALTALT 3

// Parameter thresholds for FLAGs
#define LOW_AC_THRESHOLD        5
#define LONG_RANGE_THRESHOLD    500e3
#define MINIMUM_ALLOWED_ALLELES 5		// Minimum number of alleles required for math to work out in the unphased case
#define ALLOWED_ROUNDING_ERROR  0.001

#define TWK_HAP_FREQ(A,POS) ((double)(A).alleleCounts[POS] / (A).totalHaplotypeCounts)

#if SLAVE_DEBUG_MODE == 3
static uint32_t twk_debug_pos1 = 0;
static uint32_t twk_debug_pos1_2 = 0;
static uint32_t twk_debug_pos2 = 0;
static uint32_t twk_debug_pos2_2 = 0;
#endif

/****************************
*  SIMD definitions
****************************/
#if SIMD_AVAILABLE == 1

#ifdef _mm_popcnt_u64
#define POPCOUNT_ITER	_mm_popcnt_u64
#else
#define POPCOUNT_ITER	__builtin_popcountll
#endif

#define UNPHASED_UPPER_MASK	170  // 10101010b
#define UNPHASED_LOWER_MASK	85   // 01010101b
#define FILTER_UNPHASED_BYTE(A, B)            ((((((A) & UNPHASED_UPPER_MASK) | ((B) & UNPHASED_LOWER_MASK)) & UNPHASED_LOWER_MASK) << 1) & (A))
#define FILTER_UNPHASED_BYTE_PAIR(A, B, C, D) ((FILTER_UNPHASED_BYTE((A), (B)) >> 1) | FILTER_UNPHASED_BYTE((C), (D)))
#define FILTER_UNPHASED_BYTE_SPECIAL(A)       ((((A) >> 1) & (A)) & UNPHASED_LOWER_MASK)

#if SIMD_VERSION == 6 // AVX-512: UNTESTED
#define VECTOR_TYPE	__m512i
const VECTOR_TYPE ONE_MASK         = _mm512_set1_epi8(255); // 11111111b
const VECTOR_TYPE maskUnphasedHigh = _mm512_set1_epi8(UNPHASED_UPPER_MASK);	// 10101010b
const VECTOR_TYPE maskUnphasedLow  = _mm512_set1_epi8(UNPHASED_LOWER_MASK);	// 01010101b

#define PHASED_ALTALT(A,B)        _mm512_and_si512(A, B)
#define PHASED_REFREF(A,B)        _mm512_and_si512(_mm512_xor_si512(A, ONE_MASK), _mm512_xor_si512(B, ONE_MASK))
#define PHASED_ALTREF(A,B)        _mm512_and_si512(_mm512_xor_si512(A, B), B)
#define PHASED_REFALT(A,B)        _mm512_and_si512(_mm512_xor_si512(A, B), A)
#define PHASED_ALTALT_MASK(A,B,M) _mm512_and_si512(PHASED_ALTALT(A, B), M)
#define PHASED_REFREF_MASK(A,B,M) _mm512_and_si512(PHASED_REFREF(A, B), M)
#define PHASED_ALTREF_MASK(A,B,M) _mm512_and_si512(PHASED_ALTREF(A, B), M)
#define PHASED_REFALT_MASK(A,B,M) _mm512_and_si512(PHASED_REFALT(A, B), M)
#define MASK_MERGE(A,B)           _mm512_xor_si512(_mm512_or_si512(A, B), ONE_MASK)

#define POPCOUNT(A, B) {                                \
	__m256i tempA = _mm512_extracti64x4_epi64(B, 0);    \
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 0)); \
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 1)); \
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 2)); \
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 3)); \
	tempA = _mm512_extracti64x4_epi64(B, 1);            \
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 0)); \
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 1)); \
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 2)); \
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 3)); \
}

#define FILTER_UNPHASED(A, B)            _mm512_and_si512(_mm512_slli_epi64(_mm512_and_si512(_mm512_or_si512(_mm512_and_si512(A, maskUnphasedHigh),_mm512_and_si512(B, maskUnphasedLow)), maskUnphasedLow), 1), A)
#define FILTER_UNPHASED_PAIR(A, B, C, D) _mm512_or_si512(_mm512_srli_epi64(FILTER_UNPHASED(A, B), 1), FILTER_UNPHASED(C, D))
#define FILTER_UNPHASED_SPECIAL(A)       _mm512_and_si512(_mm512_and_si512(_mm512_srli_epi64(A, 1), A), maskUnphasedLow)


#elif SIMD_VERSION == 5 // AVX2
#define VECTOR_TYPE	__m256i
const VECTOR_TYPE ONE_MASK         = _mm256_set1_epi8(255); // 11111111b
const VECTOR_TYPE maskUnphasedHigh = _mm256_set1_epi8(UNPHASED_UPPER_MASK);	// 10101010b
const VECTOR_TYPE maskUnphasedLow  = _mm256_set1_epi8(UNPHASED_LOWER_MASK);	// 01010101b

#define PHASED_ALTALT(A,B)        _mm256_and_si256(A, B)
#define PHASED_REFREF(A,B)        _mm256_and_si256(_mm256_xor_si256(A, ONE_MASK), _mm256_xor_si256(B, ONE_MASK))
#define PHASED_ALTREF(A,B)        _mm256_and_si256(_mm256_xor_si256(A, B), B)
#define PHASED_REFALT(A,B)        _mm256_and_si256(_mm256_xor_si256(A, B), A)
#define PHASED_ALTALT_MASK(A,B,M) _mm256_and_si256(PHASED_ALTALT(A, B), M)
#define PHASED_REFREF_MASK(A,B,M) _mm256_and_si256(PHASED_REFREF(A, B), M)
#define PHASED_ALTREF_MASK(A,B,M) _mm256_and_si256(PHASED_ALTREF(A, B), M)
#define PHASED_REFALT_MASK(A,B,M) _mm256_and_si256(PHASED_REFALT(A, B), M)
#define MASK_MERGE(A,B)           _mm256_xor_si256(_mm256_or_si256(A, B), ONE_MASK)

// Software intrinsic popcount
#define POPCOUNT(A, B) {							\
	A += POPCOUNT_ITER(_mm256_extract_epi64(B, 0));	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(B, 1));	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(B, 2));	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(B, 3));	\
}

#define FILTER_UNPHASED(A, B)            _mm256_and_si256(_mm256_slli_epi64(_mm256_and_si256(_mm256_or_si256(_mm256_and_si256(A, maskUnphasedHigh),_mm256_and_si256(B, maskUnphasedLow)), maskUnphasedLow), 1), A)
#define FILTER_UNPHASED_PAIR(A, B, C, D) _mm256_or_si256(_mm256_srli_epi64(FILTER_UNPHASED(A, B), 1), FILTER_UNPHASED(C, D))
#define FILTER_UNPHASED_SPECIAL(A)       _mm256_and_si256(_mm256_and_si256(_mm256_srli_epi64(A, 1), A), maskUnphasedLow)

#elif SIMD_VERSION >= 2 // SSE2+
#define VECTOR_TYPE	__m128i
const VECTOR_TYPE ONE_MASK         = _mm_set1_epi8(255); // 11111111b
const VECTOR_TYPE maskUnphasedHigh = _mm_set1_epi8(UNPHASED_UPPER_MASK);	// 10101010b
const VECTOR_TYPE maskUnphasedLow  = _mm_set1_epi8(UNPHASED_LOWER_MASK);	// 01010101b

#define PHASED_ALTALT(A,B)        _mm_and_si128(A, B)
#define PHASED_REFREF(A,B)        _mm_and_si128(_mm_xor_si128(A, ONE_MASK), _mm_xor_si128(B, ONE_MASK))
#define PHASED_ALTREF(A,B)        _mm_and_si128(_mm_xor_si128(A, B), B)
#define PHASED_REFALT(A,B)        _mm_and_si128(_mm_xor_si128(A, B), A)
#define PHASED_ALTALT_MASK(A,B,M) _mm_and_si128(PHASED_ALTALT(A, B), M)
#define PHASED_REFREF_MASK(A,B,M) _mm_and_si128(PHASED_REFREF(A, B), M)
#define PHASED_ALTREF_MASK(A,B,M) _mm_and_si128(PHASED_ALTREF(A, B), M)
#define PHASED_REFALT_MASK(A,B,M) _mm_and_si128(PHASED_REFALT(A, B), M)
#define MASK_MERGE(A,B)           _mm_xor_si128(_mm_or_si128(A, B), ONE_MASK)

__attribute__((always_inline))
static inline void popcnt128(uint64_t& a, const __m128i n) {
	a += _mm_popcnt_u64(_mm_cvtsi128_si64(n)) + _mm_popcnt_u64(_mm_cvtsi128_si64(_mm_unpackhi_epi64(n, n)));
}

#if SIMD_VERSION >= 3
#define POPCOUNT(A, B) {								\
	A += POPCOUNT_ITER(_mm_extract_epi64(B, 0));		\
	A += POPCOUNT_ITER(_mm_extract_epi64(B, 1));		\
}
#else
#define POPCOUNT(A, B) {      \
	uint64_t temp = _mm_extract_epi16(B, 0) << 6 | _mm_extract_epi16(B, 1) << 4 | _mm_extract_epi16(B, 2) << 2 | _mm_extract_epi16(B, 3); \
	A += POPCOUNT_ITER(temp); \
	temp = _mm_extract_epi16(B, 4) << 6 | _mm_extract_epi16(B, 5) << 4 | _mm_extract_epi16(B, 6) << 2 | _mm_extract_epi16(B, 7); \
	A += POPCOUNT_ITER(temp); \
}
#endif

#define FILTER_UNPHASED(A, B)            _mm_and_si128(_mm_slli_epi64(_mm_and_si128(_mm_or_si128(_mm_and_si128(A, maskUnphasedHigh),_mm_and_si128(B, maskUnphasedLow)), maskUnphasedLow), 1), A)
#define FILTER_UNPHASED_PAIR(A, B, C, D) _mm_or_si128(_mm_srli_epi64(FILTER_UNPHASED(A, B), 1), FILTER_UNPHASED(C, D))
#define FILTER_UNPHASED_SPECIAL(A)       _mm_and_si128(_mm_and_si128(_mm_srli_epi64(A, 1), A), maskUnphasedLow)
#endif
#endif // ENDIF SIMD_AVAILABLE == 1

// Supportive structure for timings
struct twk_ld_perf {
	uint64_t* cycles;
	uint64_t* freq;
};

/**<
 * Settings/parameters for both `twk_ld` and `twk_ld_engine`.
 */
struct twk_ld_settings {
	twk_ld_settings() :
		square(true), window(false), low_memory(false), bitmaps(false),
		force_phased(false), forced_unphased(false), force_cross_intervals(false),
		c_level(1), bl_size(500), b_size(10000), l_window(1000000),
		n_threads(std::thread::hardware_concurrency()), cycle_threshold(0),
		ldd_load_type(TWK_LDD_ALL), out("-"),
		minP(1), minR2(0.1), maxR2(100), minDprime(0), maxDprime(100),
		n_chunks(1), c_chunk(0)
	{}

	std::string GetString() const{
		std::string s =  "square=" + std::string((square ? "TRUE" : "FALSE"))
		              + ",window=" + std::string((window ? "TRUE" : "FALSE"))
		              + ",low_memory=" + std::string((low_memory ? "TRUE" : "FALSE"))
		              + ",bitmaps=" + std::string((bitmaps ? "TRUE" : "FALSE"))
		              + ",force_phased=" + std::string((force_phased ? "TRUE" : "FALSE"))
		              + ",force_unphased=" + std::string((forced_unphased ? "TRUE" : "FALSE"))
		              + ",compression_level=" + std::to_string(c_level)
		              + ",block_size=" + std::to_string(bl_size)
		              + ",output_block_size=" + std::to_string(b_size)
		              + (window ? std::string(",window_size=") + std::to_string(l_window) : "")
		              + ",minP=" + std::to_string(minP)
		              + ",minR2=" + std::to_string(minR2)
		              + ",maxR2=" + std::to_string(maxR2)
		              + ",minDprime=" + std::to_string(minDprime)
		              + ",maxDprime=" + std::to_string(maxDprime)
		              + ",n_chunks=" + std::to_string(n_chunks)
		              + ",c_chunk=" + std::to_string(c_chunk)
		              + ",n_threads=" + std::to_string(n_threads)
		              + ",ldd_type=" + std::to_string((int)ldd_load_type)
		              + ",cycle_threshold=" + std::to_string(cycle_threshold);
		return(s);
	}

	bool square, window, low_memory, bitmaps; // using square compute, using window compute
	bool force_phased, forced_unphased, force_cross_intervals;
	int32_t c_level, bl_size, b_size, l_window; // compression level, block_size, output block size, window size in bp
	int32_t n_threads, cycle_threshold, ldd_load_type;
	std::string in, out; // input file, output file/cout
	double minP, minR2, maxR2, minDprime, maxDprime;
	int32_t n_chunks, c_chunk;
	std::vector<std::string> ival_strings; // unparsed interval strings
};

/**<
 * Load balancer for calculating linkage-disequilibrium. Partitions the total
 * problem into psuedo-balanced subproblems. The size and number of sub-problems
 * can be parameterized.
 */
struct twk_ld_balancer {
	twk_ld_balancer() : diag(false), n(0), p(0), c(0), fromL(0), toL(0), fromR(0), toR(0), n_m(0){}

	bool BuildCrossSequence(uint32_t n_blocks,
	                        uint32_t desired_parts,
	                        uint32_t chosen_part)
	{
		// compare sections of blocks
	}

	/**<
	 * Find the desired target subproblem range as a tuple (fromL,toL,fromR,toR).
	 * @param n_blocks      Total number of blocks.
	 * @param desired_parts Desired number of subproblems to solve.
	 * @param chosen_part   Target subproblem we are interested in getting the ranges for.
	 * @return              Return TRUE upon success or FALSE otherwise.
	 */
	bool Build(uint32_t n_blocks,
	           uint32_t desired_parts,
	           uint32_t chosen_part)
	{
		if(chosen_part >= desired_parts){
			std::cerr << utility::timestamp("ERROR","BALANCER") << "Illegal chosen block: " << chosen_part << " >= " << desired_parts << std::endl;
			return false;
		}

		n = n_blocks; p = desired_parts; c = chosen_part;
		if(p > n){
			std::cerr << utility::timestamp("ERROR","BALANCER") << "Illegal desired number of blocks! You are asking for more subproblems than there are blocks available (" << p << ">" << n << ")..." << std::endl;
			return false;
		}

		if(p == 1){
			p = 1; c = 0;
			fromL = 0; toL = n_blocks; fromR = 0; toR = n_blocks;
			n_m = n_blocks; diag = true;
			return true;
		}

		uint32_t factor = 0;
		for(uint32_t i = 1; i < desired_parts; ++i){
			if( (((i*i) - i) / 2) + i == desired_parts ){
				//std::cerr << "factor is " << i << std::endl;
				factor = i;
				break;
			}
		}

		if(factor == 0){
			std::cerr << utility::timestamp("ERROR","BALANCER") << "Could not partition into " << desired_parts << " number of subproblems. This number is not a function of x!2 + x..." << std::endl;
			return false;
		}

		// cycle
		uint32_t chunk_size = n / factor;
		uint32_t fL = 0, tL = 0, fR = 0, tR = 0;
		for(uint32_t i = 0, k = 0; i < factor; ++i){ // rows
			for(uint32_t j = i; j < factor; ++j, ++k){ // cols
				tR = (j + 1 == factor ? n_blocks : chunk_size*(j+1));
				fR = tR - chunk_size;
				tL = (i + 1 == factor ? n_blocks : chunk_size*(i+1));
				fL = tL - chunk_size;

				//std::cerr << fL << "-" << tL << "->" << fR << "-" << tR << " total=" << n_blocks << " chunk=" << chunk_size << "desired=" << desired_parts << std::endl;
				if(k == chosen_part){
					//std::cerr << "chosen part:" << std::endl;
					fromL = fL; toL = tL; fromR = fR; toR = tR;
					n_m = (toL - fromL) + (toR - fromR); diag = false;
					if(i == j){ n_m = toL - fromL; diag = true; }
					return true;
				}
			}
		}
		return true;
	}

public:
	bool diag; // is selectd chunk diagonal
	uint32_t n, p, c; // number of available blocks, desired parts, chosen part
	uint32_t fromL, toL, fromR, toR;
	uint32_t n_m; // actual blocks used
};

/**<
 * Work balancer for twk_ld_engine threads. Uses a non-blocking spinlock to produce a
 * tuple (from,to) of integers representing the start ref block and dst block.
 * This approach allows perfect load-balancing at a small overall CPU cost.
 */
struct twk_ld_ticker {
	typedef bool (twk_ld_ticker::*get_func)(uint32_t& from, uint32_t& to, uint8_t& type);

	twk_ld_ticker() :
		diag(false), window(false),
		n_perf(0), i(0), j(0),
		fL(0), tL(0), fR(0), tR(0), l_window(0),
		ldd(nullptr),
		_getfunc(&twk_ld_ticker::GetBlockPair)
	{}
	~twk_ld_ticker(){}

	void operator=(const twk_ld_balancer& balancer){
		diag = balancer.diag;
		fL   = balancer.fromL;
		tL   = balancer.toL;
		fR   = balancer.fromR;
		tR   = balancer.toR;
		i    = balancer.fromL;
		j    = balancer.fromR;
	}

	void operator=(const twk_ld_settings& settings){
		window   = settings.window;
		l_window = settings.l_window;
		SetWindow(settings.window);
	}

	/**<
	 * Parameterisation of window mode: should we check if blocks overlap given
	 * some maximum distance between pairs? See function `GetBlockWindow` for
	 * additional details.
	 * @param yes Set or unset window mode.
	 */
	inline void SetWindow(const bool yes = true){
		window = yes;
		_getfunc = (window ? &twk_ld_ticker::GetBlockWindow : &twk_ld_ticker::GetBlockPair);
	}

	/**<
	 * Indirection using functional pointer to actual function used. This
	 * allows us to use a singular function without writing multiple
	 * versions of downstream functions.
	 * @param from Row position
	 * @param to   Column position
	 * @param type Diagonal (1) or square (0)
	 * @return     Returns TRUE if it is possible to retrieve a new (x,y)-pair or FALSE otherwise.
	 */
	inline bool Get(uint32_t& from, uint32_t& to, uint8_t& type){ return((this->*_getfunc)(from, to, type)); }

	/**<
	 * Retrieves (x,y)-coordinates from the selected load-balancing subproblem.
	 * This variation also checks if the two blocks (x,y) can have any overlapping
	 * regions given some parameterized maximum distance.
	 * Uses a spin-lock to make this function thread-safe.
	 * @param from Row position
	 * @param to   Column position
	 * @param type Diagonal (1) or square (0)
	 * @return     Returns TRUE if it is possible to retrieve a new (x,y)-pair or FALSE otherwise.
	 */
	bool GetBlockWindow(uint32_t& from, uint32_t& to, uint8_t& type){
		spinlock.lock();

		if(j == tR){
			++i; j = (diag ? i : fR); from = i; to = j; type = 1; ++j;
			if(i == tL){ spinlock.unlock(); return false; }
			++n_perf;
			spinlock.unlock();
			return true;
		}
		if(i == tL){ spinlock.unlock(); return false; }

		// First in tgt block - last in ref block
		if(i != j){
			// check if this (x,y) pair have any overlapping intervals.
			if(ldd[j].blk->rcds[0].pos - ldd[i].blk->rcds[ldd[i].n_rec-1].pos > l_window){
				++i; j = (diag ? i : fR); from = i; to = j; type = 1; ++j;
				spinlock.unlock();
				return true;
			}
		}

		type = (i == j); from = i; to = j;
		++j; ++n_perf;

		spinlock.unlock();

		return true;
	}

	/**<
	 * Retrieves (x,y)-coordinates from the selected load-balancing subproblem.
	 * Uses a spin-lock to make this function thread-safe.
	 * @param from Row position
	 * @param to   Column position
	 * @param type Diagonal (1) or square (0)
	 * @return     Returns TRUE if it is possible to retrieve a new (x,y)-pair or FALSE otherwise.
	 */
	bool GetBlockPair(uint32_t& from, uint32_t& to, uint8_t& type){
		spinlock.lock();

		if(j == tR){ // if current position is at the last column
			++i; j = (diag ? i : fR); from = i; to = j; type = 1; ++j;
			// if current position is at the last row
			if(i == tL){ spinlock.unlock(); return false; }
			++n_perf;
			spinlock.unlock();
			return true;
		}
		// if current position is at the last row
		if(i == tL){ spinlock.unlock(); return false; }
		type = (i == j); from = i; to = j;
		++j; ++n_perf;

		spinlock.unlock();

		return true;
	}

public:
	bool diag, window;
	uint32_t n_perf, i,j;
	uint32_t fL, tL, fR, tR, l_window;
	twk1_ldd_blk* ldd;
	get_func _getfunc;
	SpinLock spinlock;
};

/**<
 * Progress ticker for calculating linkage-disequilbirium. Spawns and detaches
 * a thread to tick occasionally in the background. Slaves computing linkage-
 * disequilibrium send their progress to this ticker that collates and summarize
 * that data.
 */
struct twk_ld_progress {
	twk_ld_progress() :
		is_ticking(false), n_s(0), n_cmps(0), n_var(0),
		n_pair(0), n_out(0), b_out(0), thread(nullptr)
	{}
	~twk_ld_progress() = default;

	/**<
	 * Starts the progress ticker. Spawns a detached thread ticking every 30 seconds
	 * in the background until the flag `is_ticking` is set to FALSE or the program
	 * finishes.
	 * @return Returns a pointer to the detached thread.
	 */
	std::thread* Start(){
		delete thread;
		is_ticking = true;
		thread = new std::thread(&twk_ld_progress::StartTicking, this);
		thread->detach();
		return(thread);
	}

	/**<
	 * Internal function displaying the progress message every 30 seconds. This
	 * function is called exclusively by the detached thread.
	 */
	void StartTicking(){
		timer.Start();

		uint64_t variant_overflow = 99E9, genotype_overflow = 999E12;
		uint8_t variant_width = 15, genotype_width = 20;

		//char support_buffer[256];
		std::cerr << utility::timestamp("PROGRESS")
				<< std::setw(12) << "Time elapsed"
				<< std::setw(variant_width) << "Variants"
				<< std::setw(genotype_width) << "Genotypes"
				<< std::setw(15) << "Output"
				<< std::setw(10) << "Progress"
				<< "\tEst. Time left" << std::endl;

		std::this_thread::sleep_for(std::chrono::seconds(30)); // first sleep
		while(is_ticking){
			// Triggered every cycle (119 ms)
			//if(this->Detailed)
			//	this->GetElapsedTime();

			//if(i % 252 == 0){ // Approximately every 30 sec (30e3 / 119 = 252)
				//const double ComparisonsPerSecond = (double)n_var.load()/timer.Elapsed().count();
				if(n_var.load() > variant_overflow){ variant_width += 3; variant_overflow *= 1e3; }
				if(n_var.load() > genotype_overflow){ genotype_width += 3; genotype_overflow *= 1e3; }


				//const uint32_t n_p = sprintf(&support_buffer[0], "%0.3f", 0);
				//support_buffer[n_p] = '%';
				if(n_cmps){
					std::cerr << utility::timestamp("PROGRESS")
							<< std::setw(12) << timer.ElapsedString()
							<< std::setw(variant_width) << utility::ToPrettyString(n_var.load())
							<< std::setw(genotype_width) << utility::ToPrettyString(n_var.load()*n_s)
							<< std::setw(15) << utility::ToPrettyString(n_out.load())
							<< std::setw(10) << (double)n_var.load()/n_cmps*100 << "%\t"
							<< (double)0 << std::endl;
				} else {
					std::cerr << utility::timestamp("PROGRESS")
							<< std::setw(12) << timer.ElapsedString()
							<< std::setw(variant_width) << utility::ToPrettyString(n_var.load())
							<< std::setw(genotype_width) << utility::ToPrettyString(n_var.load()*n_s)
							<< std::setw(15) << utility::ToPrettyString(n_out.load())
							<< std::setw(10) << 0 << '\t'
							<< 0 << std::endl;
				}
			//}
			std::this_thread::sleep_for(std::chrono::seconds(30));
		}
	}

	/**<
	 * Print out the final tally of time elapsed, number of variants computed,
	 * and average throughput. This method cannot be made const as the function
	 * ElapsedString in the Timer class internally updates a buffer for performance
	 * reasons. This has no consequence as this function is ever only called once.
	 */
	void PrintFinal(){
		std::cerr << utility::timestamp("PROGRESS") << this->timer.ElapsedString() << "\t" << utility::ToPrettyString(n_var.load()) << "\t" << utility::ToPrettyString(n_var.load()*n_s) << "\t" << utility::ToPrettyString(n_out.load()) << std::endl;
		std::cerr << utility::timestamp("PROGRESS") << utility::ToPrettyString((uint64_t)((double)n_var.load()/timer.Elapsed().count())) << "\t" << utility::ToPrettyString((uint64_t)(((double)n_var.load()*n_s)/timer.Elapsed().count())) << std::endl;
		std::cerr << utility::timestamp("PROGRESS") << "Finished" << std::endl;
	}

public:
	bool is_ticking;
	uint32_t n_s; // number of samples
	uint64_t n_cmps; // number of comparisons we estimate to perform
	std::atomic<uint64_t> n_var, n_pair, n_out, b_out; // counters used by ld threads
	std::thread* thread; // detached thread
	Timer timer; // timer instance
};


/**<
 * Supportive counter struct for when using the vectorized instructions
 * for computing linkage-disequilibrium.
 */
struct twk_ld_simd {
public:
	twk_ld_simd(void);
	~twk_ld_simd();

public:
	uint64_t *counters;
	uint8_t  *scalarA, *scalarB, *scalarC, *scalarD;
} __attribute__((aligned(16)));

/**<
 * Generic supportive counter struct for when calculating linkage-disequilibrium.
 * This struct have internal counters for haplotype and allele counts and for the
 * total number of allele/haplotype counts.
 */
struct twk_ld_count {
	twk_ld_count(): totalHaplotypeCounts(0)
	{
		// Initialize counters to 0. This is generally not necessary as each
		// function computing LD clears these. However, it is good practice.
		memset(alleleCounts, 171, sizeof(uint64_t)*171);
		memset(haplotypeCounts, 4, sizeof(uint64_t)*4);
	}
	~twk_ld_count(){}

	void resetPhased(void){
		this->alleleCounts[0]  = 0;
		this->alleleCounts[1]  = 0;
		this->alleleCounts[4]  = 0;
		this->alleleCounts[5]  = 0;
		haplotypeCounts[0] = 0;
		haplotypeCounts[1] = 0;
		haplotypeCounts[2] = 0;
		haplotypeCounts[3] = 0;
		// All other values can legally overflow
		// They are not used
	}

	void resetUnphased(void){
		this->alleleCounts[0]  = 0;
		this->alleleCounts[1]  = 0;
		this->alleleCounts[4]  = 0;
		this->alleleCounts[5]  = 0;
		this->alleleCounts[16] = 0;
		this->alleleCounts[17] = 0;
		this->alleleCounts[20] = 0;
		this->alleleCounts[21] = 0;
		this->alleleCounts[64] = 0;
		this->alleleCounts[65] = 0;
		this->alleleCounts[68] = 0;
		this->alleleCounts[69] = 0;
		this->alleleCounts[80] = 0;
		this->alleleCounts[81] = 0;
		this->alleleCounts[84] = 0;
		this->alleleCounts[85] = 0;
		// All other values can legally overflow
		// They are not used
	}

	// Counters
	uint64_t alleleCounts[171]; // equivalent to 10101010b = (MISS,MISS),(MISS,MISS)
	uint64_t haplotypeCounts[4]; // haplotype counts
	uint64_t totalHaplotypeCounts; // Total number of alleles
};

/**<
 * Internal engine for computing linkage-disequilibrium. Invocation of these
 * functions are generally performed outside of this class in the spawned
 * slaves. This separation of invokation and compute allows for embarassingly
 * parallel compute of non-overlapping blocks without blocking resources.
 */
class twk_ld_engine {
public:
	struct phased_helper {
		inline double& operator[](const int p){ return(vals[p]); }
		double vals[8];
	};

	phased_helper phased_help;

	// Function definitions to computing functions defined below.
	typedef bool (twk_ld_engine::*func)(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf);
	typedef bool (twk_ld_engine::*ep[10])(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf);

public:
	twk_ld_engine() :
		n_samples(0), n_out(0), n_lim(10000), n_out_tick(250),
		byte_width(0), byte_aligned_end(0), vector_cycles(0),
		phased_unbalanced_adjustment(0), unphased_unbalanced_adjustment(0), t_out(0),
		index(nullptr), writer(nullptr), progress(nullptr), list_out(nullptr)
	{
		memset(n_method, 0, sizeof(uint64_t)*10);
	}

	~twk_ld_engine(){ delete[] list_out; }

	/**<
	 * Set the number of samples in the target file. This function is mandatory
	 * to call prior to computing as it sets essential paramters such as the number
	 * of samples, the byte alignment of vectorized bit-vectors and the adjustment
	 * given the unused overhangs.
	 * @param samples Number of samples in the target file.
	 */
	void SetSamples(const uint32_t samples){
		n_samples  = samples;
		byte_width = ceil((double)samples/4);
		byte_aligned_end = byte_width/(GENOTYPE_TRIP_COUNT/4)*(GENOTYPE_TRIP_COUNT/4);
		vector_cycles    = byte_aligned_end*4/GENOTYPE_TRIP_COUNT;
		phased_unbalanced_adjustment   = (samples*2)%8;
		unphased_unbalanced_adjustment = samples%4;
	}

	void SetBlocksize(const uint32_t s){
		assert(s % 2 == 0);
		blk_f.clear(); blk_r.clear();
		blk_f.resize(s + 100);
		blk_r.resize(s + 100);
		obuf.resize(s * sizeof(twk1_two_t));
		ibuf.resize(s * sizeof(twk1_two_t));

		n_out = 0; n_lim = s;
	}

	/**<
	 * Phased/unphased functions for calculating linkage-disequilibrium. These
	 * functions are all prefixed with Phased_ or Unphased_. All these functions
	 * share the same interface in order to allow them being targetted by a shared
	 * function pointer.
	 * @param b1   Left twk1_ldd_blk reference.
	 * @param p1   Left relative offset into the left twk1_ldd_blk reference.
	 * @param b2   Right twk_1_ldd_blk reference.
	 * @param p2   Right relative offset into the right twk1_ldd_blk reference.
	 * @param perf Pointer to a twk_ld_perf object if performance measure are to be taken. This value can be set to nullptr.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool PhasedRunlength(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool PhasedList(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool PhasedListVector(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);

	bool PhasedListSpecial(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);

	bool PhasedVectorized(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool PhasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool PhasedBitmap(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool PhasedMath(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2);
	bool UnphasedRunlength(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool UnphasedVectorized(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool UnphasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);

	// Specialized hybrid functions
	__attribute__((always_inline))
	inline bool HybridUnphased(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr){
		if(b1.blk->rcds[p1].gt->n + b2.blk->rcds[p2].gt->n < 40)
			return(UnphasedRunlength(b1,p1,b2,p2,perf));
		else return(UnphasedVectorizedNoMissing(b1,p1,b2,p2,perf));
	}

	// Debugging function for invoking each algorithm on the same data.
	inline bool AllAlgorithms(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr){
		this->PhasedVectorized(b1,p1,b2,p2,perf);
		this->PhasedVectorizedNoMissing(b1,p1,b2,p2,perf);
		this->UnphasedVectorized(b1,p1,b2,p2,perf);
		this->UnphasedVectorizedNoMissing(b1,p1,b2,p2,perf);
		this->PhasedRunlength(b1,p1,b2,p2,perf);
		this->PhasedList(b1,p1,b2,p2,perf);
		this->UnphasedRunlength(b1,p1,b2,p2,perf);
		return true;
	}

	// Unphased math
	bool UnphasedMath(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2);
	double ChiSquaredUnphasedTable(const double target, const double p, const double q);
	bool ChooseF11Calculate(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2,const double target, const double p, const double q);

	/**<
	 * Compress the internal output buffers consisting of twk1_two_t records.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool CompressBlock();
	bool CompressFwd();
	bool CompressRev();

public:
	uint32_t n_samples;
	uint32_t n_out, n_lim, n_out_tick;
	uint32_t byte_width; // Number of bytes required per variant site
	uint32_t byte_aligned_end; // End byte position
	uint32_t vector_cycles; // Number of SIMD cycles (genotypes/2/vector width)
	uint32_t phased_unbalanced_adjustment; // Modulus remainder
	uint32_t unphased_unbalanced_adjustment; // Modulus remainder
	uint64_t t_out; // number of bytes written
	uint64_t n_method[10];

	IndexEntryOutput irecF, irecR;
	ZSTDCodec zcodec; // reusable zstd codec instance with internal context.

	twk_ld_settings settings;
	twk1_two_block_t blk_f, blk_r;
	twk_buffer_t ibuf, obuf; // buffers
	twk_ld_simd  helper_simd; // simd counters
	twk_ld_count helper; // regular counters
	twk1_two_t cur_rcd;
	twk1_two_t cur_rcd2;
	IndexOutput* index;
	twk_writer_t* writer;
	twk_ld_progress* progress;
	uint32_t* list_out;
};

struct twk_ld_slave {
	twk_ld_slave() : n_s(0), n_total(0),
		i_start(0), j_start(0), prev_i(0), prev_j(0), n_cycles(0),
		ticker(nullptr), thread(nullptr), ldd(nullptr),
		progress(nullptr), settings(nullptr), itree(nullptr)
	{}

	~twk_ld_slave(){ delete thread; }

	/**<
	 * Primary subroutine for starting the slave thread to compute its designated
	 * region(s) of linkage-disequilibrium.
	 * @return Returns a pointer to the spawned slave thread.
	 */
	std::thread* Start(){
		delete thread;

		if(settings->force_phased && settings->low_memory && settings->bitmaps)
			if(settings->window) thread = new std::thread(&twk_ld_slave::CalculatePhasedBitmapWindow, this, nullptr);
			else thread = new std::thread(&twk_ld_slave::CalculatePhasedBitmap, this, nullptr);
		else if(settings->force_phased){
			if(settings->window) thread = new std::thread(&twk_ld_slave::CalculatePhasedWindow, this, nullptr);
			else thread = new std::thread(&twk_ld_slave::CalculatePhased, this, nullptr);
		} else if(settings->forced_unphased){
			thread = new std::thread(&twk_ld_slave::CalculateUnphased, this, nullptr);
		}
		else {
			thread = new std::thread(&twk_ld_slave::Calculate, this, nullptr);
		}

		return(thread);
	}

	/**<
	 * Update local blocks. Either copies pointer to pre-allocated containers or
	 * construct new on-the-fly given the memory parametrisation. Internally
	 * ivokes the subroutines prefixed with `UpdateBlocks`.
	 * @param blks  Dst twk1_ldd_blk array of size 2.
	 * @param from  Integer start offset.
	 * @param to    Interger end offset.
	 */
	inline void UpdateBlocks(twk1_ldd_blk* blks, const uint32_t& from, const uint32_t& to){
		if(settings->low_memory == false) this->UpdateBlocksPreloaded(blks,from,to);
		else this->UpdateBlocksGenerate(blks,from,to);
	}

	/**<
	 * Subroutine for retrieving a pair of twk1_ldd_blk from the pre-loaded arrays.
	 * @param blks Dst pointers to twk1_ldd_blks.
	 * @param from Virtual offset for A.
	 * @param to   Virtual offset for B.
	 */
	void UpdateBlocksPreloaded(twk1_ldd_blk* blks, const uint32_t& from, const uint32_t& to){
		const uint32_t add = ticker->diag ? 0 : (ticker->tL - ticker->fL);
		blks[0].SetPreloaded(ldd[from - i_start]);
		blks[1].SetPreloaded(ldd[add + (to - j_start)]);
		prev_i = from - i_start; prev_j = add + (to - j_start); ++n_cycles;
	}

	/**<
	 * Subroutine for online pre-computing of a pair of twk1_ldd_blk. This
	 * subroutine is used exclusively in low-memory mode.
	 * @param blks Dst pointers to twk1_ldd_blks.
	 * @param from Virtual offset for A.
	 * @param to   Virtual offset for B.
	 */
	void UpdateBlocksGenerate(twk1_ldd_blk* blks, const uint32_t& from, const uint32_t& to){
		const uint32_t add = ticker->diag ? 0 : (ticker->tL - ticker->fL);

		if(n_cycles == 0 || prev_i != from - i_start){
			//std::cerr << "newA" << std::endl;
			//delete[] blocks[0].vec;  blocks[0].vec = nullptr;
			//delete[] blocks[0].list; blocks[0].list = nullptr;
			blks[0] = ldd[from - i_start];
			blks[0].Inflate(n_s, settings->ldd_load_type, true);
		}
		else {
			blks[0].blk = ldd[from - i_start].blk;
			blks[0].n_rec = ldd[from - i_start].blk->n;
			//std::cerr << "recycleA" << std::endl;
		}
		if(n_cycles == 0 || prev_j != add + (to - j_start)){
			//std::cerr << "newB" << std::endl;
			//delete[] blks[1].vec;  blks[1].vec = nullptr;
			//delete[] blks[1].list; blks[1].list = nullptr;
			blks[1] = ldd[add + (to - j_start)];
			blks[1].Inflate(n_s, settings->ldd_load_type, true);
		}
		else {
			blks[1].blk = ldd[add + (to - j_start)].blk;
			blks[1].n_rec = ldd[add + (to - j_start)].blk->n;
			//std::cerr << "recycleB" << std::endl;
		}

		prev_i = from - i_start; prev_j = add + (to - j_start); ++n_cycles;
	}

	/**<
	 *
	 * @param rcds0
	 * @param rcds1
	 * @param blocks
	 * @param type
	 * @param perf
	 */
	void Phased(const twk1_t* rcds0,
	            const twk1_t* rcds1,
	            const twk1_ldd_blk* blocks,
	            const uint8_t type,
	            twk_ld_perf* perf = nullptr)
	{
		// Heuristically determined linear model at varying number of samples
		// using SSE4.2
		// y = 0.008145*n_s + 32.8227
		//const uint32_t thresh_nomiss = 0.008145*n_s + 32.8227;
		// RLEP-BVP intersection
		// y = 0.0047*n_s + 5.2913
		const uint32_t thresh_miss = 0.0047*n_s + 5.2913;
		//std::cerr << thresh_nomiss << "," << thresh_miss << std::endl;

		if(type == 1){
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
					if(rcds0[i].ac + rcds0[j].ac <= 2){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						engine.PhasedListVector(blocks[0],i,blocks[0],j,perf);
						//engine.PhasedListSpecial(blocks[0],i,blocks[0],j,perf);
						/*if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss){
							engine.PhasedList(blocks[0],i,blocks[0],j,perf);
							//engine.PhasedListSpecial(blocks[0],i,blocks[0],j,perf);
							//engine.PhasedListVector(blocks[0],i,blocks[0],j,perf);
						} else {
							engine.PhasedVectorizedNoMissing(blocks[0],i,blocks[0],j,perf);
						}*/
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
						else
							engine.PhasedVectorized(blocks[0],i,blocks[0],j,perf);
					}
				}
			}
			progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
		} else {
			// Cache blocking
			const uint32_t bsize = (256e3/2) / (2*n_s/8);
			const uint32_t n_blocks1 = blocks[0].n_rec / bsize;
			const uint32_t n_blocks2 = blocks[1].n_rec / bsize;

			uint64_t d = 0;
			for(uint32_t ii = 0; ii < n_blocks1*bsize; ii += bsize){
				for(uint32_t jj = 0; jj < n_blocks2*bsize; jj += bsize){
					for(uint32_t i = ii; i < ii + bsize; ++i){
						for(uint32_t j = jj; j < jj + bsize; ++j){
							++d;
							if( rcds0[i].ac + rcds1[j].ac <= 2){
								continue;
							}

							//std::cerr << ii << "/" << n_blocks1 << "," << jj << "/" << n_blocks2 << "," << i << "/" << ii+bsize << "," << j << "/" << jj+bsize << " bsize=" << bsize << std::endl;
							if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
								engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
								//engine.PhasedListSpecial(blocks[0],i,blocks[1],j,perf);

								/*if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss){
									engine.PhasedList(blocks[0],i,blocks[1],j,perf);
									//engine.PhasedListSpecial(blocks[0],i,blocks[1],j,perf);
									//engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
								} else {
									engine.PhasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
								}*/
							} else {
								if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
									engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
								else
									engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
							}
						}
						//std::cerr << "end j" << std::endl;
					}
					//std::cerr << "end i" << std::endl;
				}
				// residual j that does not fit in a block
				for(uint32_t i = ii; i < ii + bsize; ++i){
					for(uint32_t j = n_blocks2*bsize; j < blocks[1].n_rec; ++j){
						++d;
						if( rcds0[i].ac + rcds1[j].ac <= 2){
							continue;
						}

						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
							//engine.PhasedListSpecial(blocks[0],i,blocks[1],j,perf);
							/*if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss){
								engine.PhasedList(blocks[0],i,blocks[1],j,perf);
								//engine.PhasedListSpecial(blocks[0],i,blocks[1],j,perf);
								//engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
							} else {
								engine.PhasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
							}*/
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
							else
								engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
						}

					}
				}
			}

			//std::cerr << "i=" << n_blocks1*bsize << "/" << blocks[0].n_rec << std::endl;
			//std::cerr << "j=" << n_blocks2*bsize << "/" << blocks[1].n_rec << std::endl;

			for(uint32_t i = n_blocks1*bsize; i < blocks[0].n_rec; ++i){
				for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
					++d;
					if( rcds0[i].ac + rcds1[j].ac <= 2){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
						//engine.PhasedListSpecial(blocks[0],i,blocks[1],j,perf);
						/*if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss){
							engine.PhasedList(blocks[0],i,blocks[1],j,perf);
							//engine.PhasedListSpecial(blocks[0],i,blocks[1],j,perf);
							//engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
						} else {
							engine.PhasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
						}*/
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
						else
							engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
					}

				}
			}

			//std::cerr << "done=" << d << "/" << blocks[0].n_rec * blocks[1].n_rec << std::endl;
			progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
		}
	}

	void Unphased(const twk1_t* rcds0, const twk1_t* rcds1, const twk1_ldd_blk* blocks, const uint8_t type, twk_ld_perf* perf = nullptr){
		// Heuristically determined linear model at varying number of samples
		// using SSE4.2
		// y = 0.0088*n_s + 18.972
		const uint32_t thresh_nomiss = 0.0088*n_s + 18.972;
		// RLEU-BVU intersection
		// y = 0.012*n_s + 22.3661
		const uint32_t thresh_miss = 0.012*n_s + 22.3661;

		if(type == 1){
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
					if(rcds0[i].ac + rcds0[j].ac <= 2){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
							engine.UnphasedRunlength(blocks[0],i,blocks[0],j,perf);
						else {
							engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[0],j,perf);
						}
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.UnphasedRunlength(blocks[0],i,blocks[0],j,perf);
						else
							engine.UnphasedVectorized(blocks[0],i,blocks[0],j,perf);
					}
				}
			}
			progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
		} else {
			const uint32_t bsize = (256e3/2) / (2*n_s/8);
			const uint32_t n_blocks1 = blocks[0].n_rec / bsize;
			const uint32_t n_blocks2 = blocks[1].n_rec / bsize;

			uint64_t d = 0;
			for(uint32_t ii = 0; ii < n_blocks1*bsize; ii += bsize){
				for(uint32_t jj = 0; jj < n_blocks2*bsize; jj += bsize){
					for(uint32_t i = ii; i < ii + bsize; ++i){
						for(uint32_t j = jj; j < jj + bsize; ++j){
							++d;
							if( rcds0[i].ac + rcds1[j].ac <= 2){
								continue;
							}

							//std::cerr << ii << "/" << n_blocks1 << "," << jj << "/" << n_blocks2 << "," << i << "/" << ii+bsize << "," << j << "/" << jj+bsize << " bsize=" << bsize << std::endl;
							if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
								if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
									engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
								else {
									engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
								}
							} else {
								if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
									engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
								else
									engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
							}
						}
						//std::cerr << "end j" << std::endl;
					}
					//std::cerr << "end i" << std::endl;
				}
				// residual j that does not fit in a block
				for(uint32_t i = ii; i < ii + bsize; ++i){
					for(uint32_t j = n_blocks2*bsize; j < blocks[1].n_rec; ++j){
						++d;
						if( rcds0[i].ac + rcds1[j].ac <= 2){
							continue;
						}

						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
								engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
							else {
								engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
							}
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
							else
								engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
						}

					}
				}
			}

			//std::cerr << "i=" << n_blocks1*bsize << "/" << blocks[0].n_rec << std::endl;
			//std::cerr << "j=" << n_blocks2*bsize << "/" << blocks[1].n_rec << std::endl;

			for(uint32_t i = n_blocks1*bsize; i < blocks[0].n_rec; ++i){
				for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
					++d;
					if( rcds0[i].ac + rcds1[j].ac <= 2){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
						else {
							engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
						}
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
						else
							engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
					}

				}
			}

			//std::cerr << "done=" << d << "/" << blocks[0].n_rec * blocks[1].n_rec << std::endl;
			progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
		}
	}

	bool CalculatePhased(twk_ld_perf* perf = nullptr){
		twk1_ldd_blk blocks[2];
		uint32_t from, to; uint8_t type;
		Timer timer; timer.Start();

		i_start = ticker->fL; j_start = ticker->fR;
		prev_i = 0; prev_j = 0;
		n_cycles = 0;
		const twk1_t* rcds0 = nullptr;
		const twk1_t* rcds1 = nullptr;

		while(true){
			if(!ticker->Get(from, to, type)) break;

			this->UpdateBlocks(blocks,from,to);
			rcds0 = blocks[0].blk->rcds;
			rcds1 = blocks[1].blk->rcds;

			// Compute phased math
			Phased(rcds0, rcds1, blocks, type, perf);
		}
		//std::cerr << "done" << std::endl;

		// if preloaded
		if(settings->low_memory == false){
			blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
			blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
		}

		// compress last
		if(this->engine.CompressBlock() == false)
			return false;

		return true;
	}

	bool CalculateUnphased(twk_ld_perf* perf = nullptr){
		twk1_ldd_blk blocks[2];
		uint32_t from, to; uint8_t type;
		Timer timer; timer.Start();

		i_start = ticker->fL; j_start = ticker->fR;
		prev_i = 0; prev_j = 0;
		n_cycles = 0;
		const twk1_t* rcds0 = nullptr;
		const twk1_t* rcds1 = nullptr;

		while(true){
			if(!ticker->Get(from, to, type)) break;

			this->UpdateBlocks(blocks,from,to);
			rcds0 = blocks[0].blk->rcds;
			rcds1 = blocks[1].blk->rcds;

			// Compute unphased math
			Unphased(rcds0, rcds1, blocks, type, perf);
		}
		//std::cerr << "done" << std::endl;

		// if preloaded
		if(settings->low_memory == false){
			blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
			blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
		}

		// compress last
		if(this->engine.CompressBlock() == false)
			return false;

		return true;
	}

	bool CalculatePhasedBitmap(twk_ld_perf* perf = nullptr){
		twk1_ldd_blk blocks[2];
		uint32_t from, to; uint8_t type;
		Timer timer; timer.Start();

		i_start = ticker->fL; j_start = ticker->fR;
		prev_i = 0; prev_j = 0;
		n_cycles = 0;

		const twk1_t* rcds0 = nullptr;
		const twk1_t* rcds1 = nullptr;


		while(true){
			if(!ticker->Get(from, to, type)) break;

			this->UpdateBlocks(blocks,from,to);
			rcds0 = blocks[0].blk->rcds;
			rcds1 = blocks[1].blk->rcds;

			if(type == 1){
				for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
					for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
						if(rcds0[i].ac + rcds0[j].ac <= 2){
							continue;
						}

						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedBitmap(blocks[0],i,blocks[0],j,perf);
						} else {
							engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
						}
					}
				}
				progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
			} else {
				for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
					for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
						if( rcds0[i].ac + rcds1[j].ac <= 2){
							continue;
						}

						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedBitmap(blocks[0],i,blocks[1],j,perf);
						} else {
							engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
						}
					}
				}
				progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
			}
		}

		// if preloaded
		if(settings->low_memory == false){
			blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
			blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
		}

		// compress last
		if(this->engine.CompressBlock() == false)
			return false;

		return true;
	}

	bool CalculatePhasedBitmapWindow(twk_ld_perf* perf = nullptr){
		twk1_ldd_blk blocks[2];
		uint32_t from, to; uint8_t type;
		Timer timer; timer.Start();

		i_start = ticker->fL; j_start = ticker->fR;
		prev_i = 0; prev_j = 0;
		n_cycles = 0;

		const twk1_t* rcds0 = nullptr;
		const twk1_t* rcds1 = nullptr;


		while(true){
			if(!ticker->Get(from, to, type)) break;

			this->UpdateBlocks(blocks,from,to);
			rcds0 = blocks[0].blk->rcds;
			rcds1 = blocks[1].blk->rcds;

			if(type == 1){
				for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
					for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
						if(blocks[0].blk->rcds[i].rid != blocks[0].blk->rcds[j].rid
						   && (blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
						{
							//std::cerr << "oor=" << blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[j].pos << std::endl;
							continue;
						}

						if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac <= 2){
							continue;
						}

						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedBitmap(blocks[0],i,blocks[0],j,perf);
						} else {
							engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
						}
					}
				}
				progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
			} else {
				for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
					for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
						if(blocks[0].blk->rcds[i].rid != blocks[1].blk->rcds[j].rid
						   && (blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
						{
							//std::cerr << "oor=" << blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
							continue;
						}

						if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac <= 2 ){
							continue;
						}

						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedBitmap(blocks[0],i,blocks[1],j,perf);
						} else {
							engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
						}
					}
				}
				progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
			}
		}

		// if preloaded
		if(settings->low_memory == false){
			blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
			blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
		}

		// compress last
		if(this->engine.CompressBlock() == false)
			return false;

		return true;
	}

	bool CalculatePhasedWindow(twk_ld_perf* perf = nullptr){
		twk1_ldd_blk blocks[2];
		uint32_t from, to; uint8_t type;
		Timer timer; timer.Start();

		i_start = ticker->fL; j_start = ticker->fR;
		prev_i = 0; prev_j = 0;
		n_cycles = 0;

		const uint32_t thresh_miss = 0.0047*n_s + 5.2913;

		const twk1_t* rcds0 = nullptr;
		const twk1_t* rcds1 = nullptr;

		while(true){
			if(!ticker->Get(from, to, type)) break;

			this->UpdateBlocks(blocks,from,to);
			rcds0 = blocks[0].blk->rcds;
			rcds1 = blocks[1].blk->rcds;

			if(type == 1){
				for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
					for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
						//std::cerr << "here=" << i << "/" << j << " diff=" << blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
 						if(blocks[0].blk->rcds[i].rid == blocks[0].blk->rcds[j].rid
						   && (blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
						{
							//std::cerr << "oor=" << blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
							//progress->n_var += i*j;
							goto end_cycle;
							//continue;
						}

						if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac <= 2){
							continue;
						}

						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedListVector(blocks[0],i,blocks[0],j,perf);
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
							else
								engine.PhasedVectorized(blocks[0],i,blocks[0],j,perf);
						}
					}
				}
				progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
			} else {
				for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
					for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
						if(blocks[0].blk->rcds[i].rid == blocks[1].blk->rcds[j].rid
						   && (blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
						{
							//std::cerr << "oor=" << blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
							goto end_cycle;
							//continue;
						}

						if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac <= 2 ){
							continue;
						}

						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
							else
								engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
						}
					}
				}
				progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
			}
			end_cycle:
			progress->n_var += 0;
		}

		// if preloaded
		if(settings->low_memory == false){
			blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
			blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
		}

		// compress last
		if(this->engine.CompressBlock() == false)
			return false;

		return true;
	}

	bool Calculate(twk_ld_perf* perf = nullptr){
		twk1_ldd_blk blocks[2];
		uint32_t from, to; uint8_t type;
		Timer timer; timer.Start();

		i_start = ticker->fL; j_start = ticker->fR;
		prev_i = 0; prev_j = 0;
		n_cycles = 0;

		// Heuristically determined
		uint32_t cycle_thresh_p = n_s / 45 == 0 ? 2 : n_s / 45;
		if(settings->cycle_threshold != 0)
			cycle_thresh_p = settings->cycle_threshold;
		uint32_t cycle_thresh_u = (n_s / 60 == 0 ? 2 : n_s / 60);
		if(settings->cycle_threshold != 0)
			cycle_thresh_u = settings->cycle_threshold;

		const uint32_t thresh_miss = 0.0047*n_s + 5.2913;

		const twk1_t* rcds0 = nullptr;
		const twk1_t* rcds1 = nullptr;

		while(true){
			if(!ticker->Get(from, to, type)) break;

			this->UpdateBlocks(blocks,from,to);
			rcds0 = blocks[0].blk->rcds;
			rcds1 = blocks[1].blk->rcds;

			if(type == 1){
				for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
					for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
						if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac <= 2){
							continue;
						}

						if(blocks[0].blk->rcds[i].an || blocks[0].blk->rcds[j].an){
							if(blocks[0].blk->rcds[i].gt->n + blocks[0].blk->rcds[j].gt->n < cycle_thresh_u){
								engine.UnphasedRunlength(blocks[0],i,blocks[0],j,nullptr);
							} else {
								engine.UnphasedVectorized(blocks[0],i,blocks[0],j,nullptr);
							}
						} else {
							if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
								engine.PhasedListVector(blocks[0],i,blocks[0],j,perf);
							} else {
								if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
									engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
								else
									engine.PhasedVectorized(blocks[0],i,blocks[0],j,perf);
							}
						}
					}
				}
				progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
			} else {
				for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
					for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
						if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac <= 2 ){
							continue;
						}

						if(blocks[0].blk->rcds[i].an || blocks[1].blk->rcds[j].an){
							if(blocks[0].blk->rcds[i].gt->n + blocks[1].blk->rcds[j].gt->n < cycle_thresh_u)
								engine.UnphasedRunlength(blocks[0],i,blocks[1],j,nullptr);
							else {
								engine.UnphasedVectorized(blocks[0],i,blocks[1],j,nullptr);
							}
						} else {
							if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
								engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
							} else {
								if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
									engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
								else
									engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
							}
						}

					}
				}
				progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
			}

		}

		// if preloaded
		if(settings->low_memory == false){
			blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
			blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
		}

		// compress last
		if(this->engine.CompressBlock() == false)
			return false;

		return true;
	}

	bool CalculatePerformance(twk_ld_engine::func f, twk_ld_perf* perf = nullptr){
		twk1_ldd_blk blocks[2];
		uint32_t from, to; uint8_t type;
		Timer timer; timer.Start();

		i_start = ticker->fL; j_start = ticker->fR;
		prev_i = 0; prev_j = 0;
		n_cycles = 0;
		const twk1_t* rcds0 = nullptr;
		const twk1_t* rcds1 = nullptr;

		while(true){
			if(!ticker->Get(from, to, type)) break;

			this->UpdateBlocks(blocks,from,to);
			rcds0 = blocks[0].blk->rcds;
			rcds1 = blocks[1].blk->rcds;

			if(type == 1){
				for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
					for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
						(engine.*f)(blocks[0],i,blocks[0],j,perf);
					}
				}
				progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
			} else {
				uint32_t bsize = (256e3/2) / (2*n_s/8);
				bsize = bsize < 1 ? 1 : bsize; // minimum block size
				const uint32_t n_blocks1 = blocks[0].n_rec / bsize;
				const uint32_t n_blocks2 = blocks[1].n_rec / bsize;

				// Cache blocking
				for(uint32_t ii = 0; ii < n_blocks1*bsize; ii += bsize){
					for(uint32_t jj = 0; jj < n_blocks2*bsize; jj += bsize){
						for(uint32_t i = ii; i < ii + bsize; ++i){
							for(uint32_t j = jj; j < jj + bsize; ++j){
								(engine.*f)(blocks[0],i,blocks[1],j,perf);
							}
						}
					}
					// residual j that does not fit in a block
					for(uint32_t i = ii; i < ii + bsize; ++i){
						for(uint32_t j = n_blocks2*bsize; j < blocks[1].n_rec; ++j){
							(engine.*f)(blocks[0],i,blocks[1],j,perf);
						}
					}
				}

				// Residual
				for(uint32_t i = n_blocks1*bsize; i < blocks[0].n_rec; ++i){
					for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
						(engine.*f)(blocks[0],i,blocks[1],j,perf);
					}
				}
				progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
			}

		}

		// if preloaded
		if(settings->low_memory == false){
			blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
			blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
		}

		return true;
	}

public:
	uint32_t n_s, n_total;
	uint32_t i_start, j_start, prev_i, prev_j, n_cycles;

	twk_ld_ticker* ticker;
	std::thread* thread;
	twk1_ldd_blk* ldd;
	twk_ld_progress* progress;
	twk_ld_settings* settings;
	twk_ld_engine engine;
	algorithm::IntervalTree<uint32_t, uint32_t>* itree;
};

// Parallel unpacker.
struct twk_ld_unpacker {
	std::thread* Start(const bool diag, const twk_ld_settings& settings){
		load = settings.ldd_load_type;
		this->diag = diag;

		bit.stream = new std::ifstream(settings.in, std::ios::binary | std::ios::in);
		if(bit.stream->good() == false){
			std::cerr << "failed to open=" << settings.in << std::endl;
			return nullptr;
		}

		if(diag) thread = new std::thread(&twk_ld_unpacker::UnpackDiagonal, this);
		else thread = new std::thread(&twk_ld_unpacker::UnpackSquare, this);
		return(thread);
	}

	bool UnpackDiagonal(){
		bit.stream->seekg(rdr->index.ent[fL].foff);
		if(bit.stream->good() == false){
			std::cerr << "failed to seek to index offset " << fL << " -> " << rdr->index.ent[fL].foff << std::endl;
			return false;
		}

		//std::cerr << "in unpack diag (" << fL << "-" << tL << ")" << std::endl;
		//std::cerr << "left-shift=" << lshift << std::endl;
		for(int i = lshift; i < lshift + (tL-fL); ++i){
			//std::cerr << "at1=" << i << "/" << lshift + (tL-fL) << std::endl;
			if(bit.NextBlock() == false){
				std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
				return false;
			}

			ldd2[i] = std::move(bit.blk);
			ldd[i].SetOwn(ldd2[i], rdr->hdr.GetNumberSamples());
			ldd[i].Inflate(rdr->hdr.GetNumberSamples(),load, resize);
		}

		std::ifstream* s = reinterpret_cast<std::ifstream*>(bit.stream);
		s->close();
		bit.stream = nullptr;

		return true;
	}

	bool UnpackSquare(){
		bit.stream->seekg(rdr->index.ent[fL].foff);
		if(bit.stream->good() == false){
			std::cerr << "failed to seek to index offset " << fL << " -> " << rdr->index.ent[fL].foff << std::endl;
			return false;
		}

		//std::cerr << "in unpack square (" << fL << "-" << tL << ")(" << fR << "-" << tR << ")" << std::endl;
		//std::cerr << "left-shift=" << lshift << std::endl;
		for(int i = lshift; i < lshift + (tL-fL); ++i){
			//std::cerr << "at1=" << i << "/" << lshift + (tL-fL) << std::endl;
			if(bit.NextBlock() == false){
				std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
				return false;
			}

			ldd2[i] = std::move(bit.blk);
			ldd[i].SetOwn(ldd2[i], rdr->hdr.GetNumberSamples());
			ldd[i].Inflate(rdr->hdr.GetNumberSamples(),load, resize);
		}

		bit.stream->seekg(rdr->index.ent[fR].foff); // seek absolute offset
		if(bit.stream->good() == false){
			std::cerr << "failed to seek to index offset " << fR << " -> " << rdr->index.ent[fR].foff << std::endl;
			return false;
		}

		// compute with local offset
		//std::cerr << "left offset=" << loff << " and right=" << roff << std::endl;
		//std::cerr << "offset=" << fR - loff << " and steps=" << loff + (tR - fR) << std::endl;
		for(int i = loff + roff; i < loff + roff + (tR - fR); ++i){
			//std::cerr << "at2=" << i << "/" << loff + roff + (tR - fR) << " with range=" << (tR - fR) << std::endl;
			if(bit.NextBlock() == false){
				std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
				return false;
			}

			ldd2[i] = std::move(bit.blk);
			ldd[i].SetOwn(ldd2[i], rdr->hdr.GetNumberSamples());
			ldd[i].Inflate(rdr->hdr.GetNumberSamples(),load, resize);
		}

		std::ifstream* s = reinterpret_cast<std::ifstream*>(bit.stream);
		s->close();
		bit.stream = nullptr;

		return true;
	}

	bool resize, diag;
	uint8_t load;
	uint32_t fL, tL, fR, tR, loff, roff, lshift;
	twk_reader* rdr;
	std::thread* thread;
	twk1_ldd_blk* ldd;
	twk1_block_t* ldd2;
	twk1_blk_iterator bit;
};

class twk_ld {
public:
	twk_ld() : n_blks(0), m_blks(0), n_vnts(0), n_tree(0), ldd(nullptr), ldd2(nullptr)
	{
	}

	~twk_ld(){
		delete[] ldd; delete[] ldd2;;
	}

	/**<
	 * Wrapper function for parsing interval strings.
	 * @param reader Reference twk_reader instance.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseIntervalStrings(twk_reader& reader){ return(intervals.ParseIntervalStrings(settings.ival_strings, reader.hdr)); }

	void operator=(const twk_ld_settings& settings){ this->settings = settings; }

	// move out
	bool WriteHeader(twk_writer_t* writer, twk_reader& reader) const{
		if(writer == nullptr)
			return false;

		tomahawk::ZSTDCodec zcodec;
		writer->write(tomahawk::TOMAHAWK_LD_MAGIC_HEADER.data(), tomahawk::TOMAHAWK_LD_MAGIC_HEADER_LENGTH);
		tomahawk::twk_buffer_t buf(256000), obuf(256000);
		buf << reader.hdr;
		if(zcodec.Compress(buf, obuf, settings.c_level) == false){
			std::cerr << "failed to compress" << std::endl;
			return false;
		}

		writer->write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		writer->write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		writer->write(obuf.data(),obuf.size());
		writer->stream.flush();
		buf.reset();

		return(writer->good());
	}

	// write out
	bool WriteFinal(twk_writer_t* writer, IndexOutput& index) const {
		if(writer == nullptr) return false;

		tomahawk::twk_buffer_t buf(256000), obuf(256000);
		tomahawk::ZSTDCodec zcodec;

		buf << index;
		if(zcodec.Compress(buf, obuf, settings.c_level) == false){
			std::cerr << "failed compression" << std::endl;
			return false;
		}

		const uint64_t offset_start_index = writer->stream.tellp();
		uint8_t marker = 0;
		writer->write(reinterpret_cast<const char*>(&marker),sizeof(uint8_t));
		writer->write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		writer->write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		writer->write(obuf.data(),obuf.size());
		writer->write(reinterpret_cast<const char*>(&offset_start_index),sizeof(uint64_t));
		writer->write(tomahawk::TOMAHAWK_FILE_EOF.data(), tomahawk::TOMAHAWK_FILE_EOF_LENGTH);
		writer->stream.flush();
		return(writer->good());
	}

	/**<
	 * Reads the desired tomahawk blocks given the balancer intervals. Internally
	 * decides if the block slicing is based on the universal set of blocks or a
	 * targetted subset (as decided by the interval slicing operation).
	 * @param reader
	 * @param bit
	 * @param balancer
	 * @param load
	 * @return
	 */
	bool LoadBlocks(twk_reader& reader,
			twk1_blk_iterator& bit,
			const twk_ld_balancer& balancer,
			const twk_ld_settings& settings)
	{
		if(intervals.overlap_blocks.size()) return(this->LoadTargetBlocks(reader, bit, balancer, settings));
		return(this->LoadAllBlocks(reader, bit, balancer, settings));
	}

	/**<
	 * Loading twk blocks for a single variant and its surrounding variants within
	 * some distance and on the same chromosome. This function loads the target variant
	 * in a single block and the other variants in separate blocks as usual. The
	 * identity of the target site is parameterized in the settings object.
	 * @param reader
	 * @param bit
	 * @param balancer
	 * @param load
	 * @return
	 */
	bool LoadTargetSingle(twk_reader& reader,
			twk1_blk_iterator& bit,
			const twk_ld_balancer& balancer,
			const uint8_t load)
	{
		if(settings.ival_strings.size() == 0){ // if have interval strings
			return false;
		}

		this->intervals.ivecs.resize(reader.hdr.GetNumberContigs());
		if(this->ParseIntervalStrings(reader) == false)
			return false;

		if(this->intervals.Build(reader.hdr.GetNumberContigs(), reader.index) == false){
			return false;
		}

		// First block has only 1 variant: the reference.
		// All other blocks
		return true;
	}

	/**<
	 * Construct interval container and trees given the pre-provided interval
	 * strings.
	 * @param reader
	 * @param bit
	 * @param load
	 * @return
	 */
	bool BuildIntervals(twk_reader& reader, twk1_blk_iterator& bit)
	{
		if(settings.ival_strings.size() == 0){ // if have interval strings
			return false;
		}

		this->intervals.ivecs.resize(reader.hdr.GetNumberContigs());
		// Parse the literal interval strings into actual interval tuples.
		if(this->ParseIntervalStrings(reader) == false)
			return false;

		// Construct interval trees and find overlaps.
		if(this->intervals.Build(reader.hdr.GetNumberContigs(), reader.index) == false){
			return false;
		}

		// Number of blocks available to iterate over.
		this->n_blks = this->intervals.overlap_blocks.size();

		return true;
	}

	/**<
	 * Loads only the target blocks that overlap with the given vector of interval
	 * tuples as parameterized in the settings object.
	 * @param reader   Reference to twk reader.
	 * @param bit      Reference to a twk block iterator.
	 * @param balancer Reference of a pre-computed load balancer.
	 * @param settings Reference of a user-paramterized settings object.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool LoadTargetBlocks(twk_reader& reader,
			twk1_blk_iterator& bit,
			const twk_ld_balancer& balancer,
			const twk_ld_settings& settings)
	{
		if(intervals.overlap_blocks.size() == 0){ // if have interval strings
			return false;
		}

		// src1 and src2 blocks are the same: e.g. (5,5) or (10,10).
		if(balancer.diag){
			//std::cerr << "is diag" << std::endl;
			//std::cerr << "load range=" << balancer.toR - balancer.fromR << std::endl;

			n_blks = balancer.toR - balancer.fromR;
			m_blks = balancer.toR - balancer.fromR;
			ldd  = new twk1_ldd_blk[m_blks];
			ldd2 = new twk1_block_t[m_blks];
		}
		// Computation is square: e.g. (5,6) or (10,21).
		else {
			//std::cerr << "is square" << std::endl;
			//std::cerr << "load range=" << (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR) << std::endl;

			n_blks = (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR);
			m_blks = (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR);
			ldd  = new twk1_ldd_blk[m_blks];
			ldd2 = new twk1_block_t[m_blks];
		}

		if(settings.low_memory)
			std::cerr << utility::timestamp("LOG") << "Running in restriced memory mode..." << std::endl;
		else
			std::cerr << utility::timestamp("LOG") << "Running in standard mode. Pre-computing data..." << std::endl;

	#if SIMD_AVAILABLE == 1
			std::cerr << utility::timestamp("LOG","SIMD") << "Vectorized instructions available: " << TWK_SIMD_MAPPING[SIMD_VERSION] << "..." << std::endl;
	#else
			std::cerr << utility::timestamp("LOG","SIMD") << "No vectorized instructions available..." << std::endl;
	#endif
			std::cerr << utility::timestamp("LOG") << "Constructing list, vector, RLE... ";

			Timer timer; timer.Start();
			if(balancer.diag){
				bit.stream->seekg(intervals.overlap_blocks[balancer.fromL]->foff);
				for(int i = 0; i < (balancer.toL - balancer.fromL); ++i){
					if(bit.NextBlock() == false){
						std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
						return false;
					}

					ldd2[i] = std::move(bit.blk);
					ldd[i].SetOwn(ldd2[i], reader.hdr.GetNumberSamples());
					ldd[i].Inflate(reader.hdr.GetNumberSamples(), settings.ldd_load_type, true);
				}
			} else {
				uint32_t offset = 0;
				bit.stream->seekg(intervals.overlap_blocks[balancer.fromL]->foff);
				for(int i = 0; i < (balancer.toL - balancer.fromL); ++i){
					if(bit.NextBlock() == false){
						std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
						return false;
					}

					ldd2[offset] = std::move(bit.blk);
					ldd[offset].SetOwn(ldd2[offset], reader.hdr.GetNumberSamples());
					ldd[offset].Inflate(reader.hdr.GetNumberSamples(),settings.ldd_load_type, true);
					++offset;
				}

				bit.stream->seekg(intervals.overlap_blocks[balancer.fromR]->foff);
				for(int i = 0; i < (balancer.toR - balancer.fromR); ++i){
					if(bit.NextBlock() == false){
						std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
						return false;
					}

					ldd2[offset] = std::move(bit.blk);
					ldd[offset].SetOwn(ldd2[offset], reader.hdr.GetNumberSamples());
					ldd[offset].Inflate(reader.hdr.GetNumberSamples(),settings.ldd_load_type, true);
					++offset;
				}
			}

			std::cerr << "Done! " << timer.ElapsedString() << std::endl;
			//std::cerr << "ldd2=" << n_blks << "/" << m_blks << std::endl;
			return(true);

		return true;
	}

	/**<
	 * Loads all available twk blocks into memory. Internally spawns the maximum
	 * possible number of unpacking threads possible (as parameterized by settings).
	 * @param reader   Reference to twk reader.
	 * @param bit      Reference to a twk block iterator.
	 * @param balancer Reference of a pre-computed load balancer.
	 * @param settings Reference of a user-paramterized settings object.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool LoadAllBlocks(twk_reader& reader,
			twk1_blk_iterator& bit,
			const twk_ld_balancer& balancer,
			const twk_ld_settings& settings)
	{
		if(reader.index.n == 0)
			return false;

		// Delete old
		delete[] ldd; delete[] ldd2;
		ldd = nullptr; ldd2 = nullptr;

		if(balancer.diag){
			//std::cerr << "is diag" << std::endl;
			//std::cerr << "load range=" << balancer.toR - balancer.fromR << std::endl;

			n_blks = balancer.toR - balancer.fromR;
			m_blks = balancer.toR - balancer.fromR;
			std::cerr << utility::timestamp("LOG") << "Allocating " << m_blks << " blocks..." << std::endl;

			ldd  = new twk1_ldd_blk[m_blks];
			ldd2 = new twk1_block_t[m_blks];
		} else {
			//std::cerr << "is square" << std::endl;
			//std::cerr << "load range=" << (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR) << std::endl;

			n_blks = (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR);
			m_blks = (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR);
			std::cerr << utility::timestamp("LOG") << "Allocating " << m_blks << " blocks..." << std::endl;

			ldd  = new twk1_ldd_blk[m_blks];
			ldd2 = new twk1_block_t[m_blks];
		}

		if(settings.low_memory)
			std::cerr << utility::timestamp("LOG") << "Running in restriced memory mode..." << std::endl;
		else
			std::cerr << utility::timestamp("LOG") << "Running in standard mode. Pre-computing data..." << std::endl;

#if SIMD_AVAILABLE == 1
		std::cerr << utility::timestamp("LOG","SIMD") << "Vectorized instructions available: " << TWK_SIMD_MAPPING[SIMD_VERSION] << "..." << std::endl;
#else
		std::cerr << utility::timestamp("LOG","SIMD") << "No vectorized instructions available..." << std::endl;
#endif

		// Distribute unpacking across multiple slaves.
		// This has insignificant performance impact on small files and/or
		// files with small number of samples (n<10,000). It has considerable
		// impact of large files when hundreds to thousands of subproblems
		// is solved.
		const uint32_t rangeL = balancer.toL - balancer.fromL;
		const uint32_t rangeR = balancer.toR - balancer.fromR;
		uint32_t unpack_threads = settings.n_threads;
		uint32_t ppthreadL = std::ceil((float)(balancer.toL - balancer.fromL) / unpack_threads);
		// Assert blocks are balanced.
		assert(balancer.toL - balancer.fromL == balancer.toR - balancer.fromR);

		if(ppthreadL == 0){
			unpack_threads = rangeL;
			//std::cerr << "changing unpack threads=" << settings.n_threads << "->" << unpack_threads << std::endl;
		} else {
			if(ppthreadL * unpack_threads > rangeL){
				unpack_threads = rangeL / ppthreadL;
				//std::cerr << "changing unpack threads=" << settings.n_threads << "->" << unpack_threads << std::endl;
			}
		}

		//std::cerr << "balance=" << (balancer.toL - balancer.fromL) << " and " << (balancer.toL - balancer.fromL) / unpack_threads << " -> " << ppthreadL << std::endl;
		std::cerr << utility::timestamp("LOG") << "Constructing list, vector, RLE..." << std::endl;
		std::cerr << utility::timestamp("LOG","THREAD") << "Unpacking using " << unpack_threads << " threads... ";
		twk_ld_unpacker* slaves = new twk_ld_unpacker[unpack_threads];
		Timer timer; timer.Start();
		for(uint32_t i = 0; i < unpack_threads; ++i){
			slaves[i].ldd = ldd;
			slaves[i].ldd2 = ldd2;
			slaves[i].rdr = &reader;
			slaves[i].resize = true;
			slaves[i].fL = balancer.fromL + ppthreadL*i; // from-left
			slaves[i].tL = i+1 == unpack_threads ? balancer.fromL+rangeL : balancer.fromL+(ppthreadL*(i+1)); // to-left
			slaves[i].fR = balancer.fromR + ppthreadL*i; // from-right
			slaves[i].tR = i+1 == unpack_threads ? balancer.fromR+rangeR : balancer.fromR+(ppthreadL*(i+1)); // to-right
			slaves[i].loff = balancer.toL - balancer.fromL; // offset from left
			slaves[i].lshift = slaves[i].fL - balancer.fromL; // offset from right
			slaves[i].roff = slaves[i].fR - balancer.fromR;
			//std::cerr << "range=" << slaves[i].fL << "->" << slaves[i].tL << " and " << slaves[i].fR << "->" << slaves[i].tR << std::endl;
		}

		for(uint32_t i = 0; i < unpack_threads; ++i) slaves[i].Start(balancer.diag, settings);
		for(uint32_t i = 0; i < unpack_threads; ++i) slaves[i].thread->join();
		delete[] slaves;

		std::cerr << "Done! " << timer.ElapsedString() << std::endl;
		return(true);
	}

	/**<
	 * Helper function to call Compute subroutine when passing a new settings
	 * object to be used.
	 * @param settings Src settings objects.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool Compute(const twk_ld_settings& settings){
		this->settings = settings;
		return(Compute());
	}

	/**<
	 * Main subroutine for computing linkage-disequilibrium as contextually
	 * determined given the user-defined parameters.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool Compute(){
		//return(ComputePerformance());

		if(settings.in.size() == 0){
			std::cerr << utility::timestamp("ERROR") << "No file-name provided..." << std::endl;
			return false;
		}

		if(settings.window && settings.n_chunks != 1){
			std::cerr << utility::timestamp("ERROR") << "Cannot use chunking in window mode!" << std::endl;
			return false;
		}

		std::cerr << utility::timestamp("LOG","READER") << "Opening " << settings.in << "..." << std::endl;

		twk_reader reader;
		if(reader.Open(settings.in) == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to open file: " << settings.in << "..." << std::endl;
			return false;
		}

		std::cerr << utility::timestamp("LOG") << "Samples: " << utility::ToPrettyString(reader.hdr.GetNumberSamples()) << "..." << std::endl;

		twk1_blk_iterator bit;
		bit.stream = reader.stream;

		Timer timer;

		//settings.ldd_load_type = TWK_LDD_ALL;
		if(settings.low_memory && settings.force_phased && settings.bitmaps){
			settings.ldd_load_type = TWK_LDD_BITMAP;
		} else if(settings.low_memory && settings.force_phased){
			settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
		} else if(settings.force_phased){
			settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
		} else if(settings.forced_unphased){
			settings.ldd_load_type = TWK_LDD_VEC;
		} else {
			settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
		}

		if(settings.ival_strings.size() == 0) n_blks = reader.index.n;
		else {
			if(this->BuildIntervals(reader, bit) == false)
				return false;
		}

		if(settings.window) settings.c_chunk = 0;

		twk_ld_balancer balancer;
		if(balancer.Build(this->n_blks, settings.n_chunks, settings.c_chunk) == false){
			return false;
		}

		std::cerr << utility::timestamp("LOG","BALANCING") << "Using ranges [" << balancer.fromL << "-" << balancer.toL << "," << balancer.fromR << "-" << balancer.toR << "] in " << (settings.window ? "window mode" : "square mode") <<"..." << std::endl;

		if(this->LoadBlocks(reader, bit, balancer, settings) == false){
			return false;
		}

		if(this->n_blks == 0){
			std::cerr << utility::timestamp("ERROR") << "No valid data available..." << std::endl;
			return true;
		}

		uint32_t n_variants = 0;
		if(balancer.diag){
			for(int i = balancer.fromL; i < balancer.toL; ++i) n_variants += reader.index.ent[i].n;
		} else {
			for(int i = balancer.fromL; i < balancer.toL; ++i) n_variants += reader.index.ent[i].n;
			for(int i = balancer.fromR; i < balancer.toR; ++i) n_variants += reader.index.ent[i].n;
		}

		std::cerr << utility::timestamp("LOG") << utility::ToPrettyString(n_variants) << " variants from " << utility::ToPrettyString(balancer.n_m) << " blocks..." << std::endl;
		std::cerr << utility::timestamp("LOG","PARAMS") << settings.GetString() << std::endl;
		std::cerr << utility::timestamp("LOG") << "Performing: " << utility::ToPrettyString(((uint64_t)n_variants * n_variants - n_variants) / 2) << " variant comparisons..." << std::endl;

		twk_ld_ticker ticker;
		ticker = balancer;
		ticker = settings;
		ticker.ldd  = ldd;

		twk_ld_progress progress;
		progress.n_s = reader.hdr.GetNumberSamples();
		if(settings.window == false){
			if(balancer.diag)
				progress.n_cmps = ((uint64_t)n_variants * n_variants - n_variants) / 2;
			else {

			}
		}
		twk_ld_slave* slaves = new twk_ld_slave[settings.n_threads];
		std::vector<std::thread*> threads(settings.n_threads);

		// Start writing file.
		twk_writer_t* writer = nullptr;
		if(settings.out.size() == 0 || (settings.out.size() == 1 && settings.out[0] == '-')){
			std::cerr << utility::timestamp("LOG","WRITER") << "Writing to " << "stdout..." << std::endl;
			writer = new twk_writer_stream;
		} else {
			std::string base_path = twk_writer_t::GetBasePath(settings.out);
			std::string base_name = twk_writer_t::GetBaseName(settings.out);
			std::string extension = twk_writer_t::GetExtension(settings.out);
			if(extension.length() == 3){
				if(strncasecmp(&extension[0], "two", 3) != 0){
					settings.out =  (base_path.size() ? base_path + "/" : "") + base_name + ".two";
				}
			} else {
				 settings.out = (base_path.size() ? base_path + "/" : "") + base_name + ".two";
			}

			std::cerr << utility::timestamp("LOG","WRITER") << "Opening " << settings.out << "..." << std::endl;
			writer = new twk_writer_file;
			if(writer->Open(settings.out) == false){
				std::cerr << utility::timestamp("ERROR", "WRITER") << "Failed to open file: " << settings.out << "..." << std::endl;
				delete writer;
				return false;
			}
		}

		// Append literal string.
		std::string calc_string = "\n##tomahawk_calcVersion=" + std::to_string(VERSION) + "\n";
		calc_string += "##tomahawk_calcCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + utility::datetime() + "\n";
		reader.hdr.literals_ += calc_string;

		if(this->WriteHeader(writer, reader) == false){
			std::cerr << utility::timestamp("ERROR","WRITER") << "Failed to write header!" << std::endl;
			return false;
		}
		// end start write

		// New index
		IndexOutput index(reader.hdr.GetNumberContigs());

		timer.Start();
		std::cerr << utility::timestamp("LOG","THREAD") << "Spawning " << settings.n_threads << " threads: ";
		for(int i = 0; i < settings.n_threads; ++i){
			slaves[i].ldd    = ldd;
			slaves[i].n_s    = reader.hdr.GetNumberSamples();
			slaves[i].ticker = &ticker;
			slaves[i].engine.SetSamples(reader.hdr.GetNumberSamples());
			slaves[i].engine.SetBlocksize(settings.b_size);
			slaves[i].engine.progress = &progress;
			slaves[i].engine.writer   = writer;
			slaves[i].engine.index    = &index;
			slaves[i].engine.settings = settings;
			slaves[i].progress = &progress;
			slaves[i].settings = &settings;
			threads[i] = slaves[i].Start();
			std::cerr << ".";
		}
		std::cerr << std::endl;

		progress.Start();

		for(int i = 0; i < settings.n_threads; ++i) threads[i]->join();
		for(int i = 0; i < settings.n_threads; ++i) slaves[i].engine.CompressBlock();
		progress.is_ticking = false;
		progress.PrintFinal();
		writer->stream.flush();

		std::cerr << utility::timestamp("LOG","THREAD") << "Thread\tOutput\tTWK-LIST\tTWK-BVP-BM\tTWK-BVP\tTWK-BVP-NM\tTWK-BVU\tTWK-BVU-NM\tTWK-RLEP\tTWK-RLEU\n";
		for(int i = 0; i < settings.n_threads; ++i){
			std::cerr << i << "\t" << utility::ToPrettyString(slaves[i].engine.n_out);
			for(int j = 0; j < 8; ++j){
				std::cerr << "\t" << utility::ToPrettyString(slaves[i].engine.n_method[j]);
			}
			std::cerr << std::endl;
		}

		//std::cerr << "performed=" << ticker.n_perf << std::endl;
		if(this->WriteFinal(writer, index) == false){
			std::cerr << utility::timestamp("ERROR","WRITER") << "Failed to write final block!" << std::endl;
			return false;
		}

		delete[] slaves; delete writer;
		std::cerr << utility::timestamp("LOG","PROGRESS") << "All done..." << timer.ElapsedString() << "!" << std::endl;

		return true;
	}

	bool ComputePerformance(){
		if(SLAVE_DEBUG_MODE != 1){
			std::cerr << utility::timestamp("ERROR","COMPILATION") << "Cannot run performance mode without compiling SLAVE_DEBUG_MODE set to 1!" << std::endl;
			return false;
		}

		if(settings.in.size() == 0){
			std::cerr << utility::timestamp("ERROR") << "No file-name provided..." << std::endl;
			return false;
		}

		if(settings.window && settings.n_chunks != 1){
			std::cerr << utility::timestamp("ERROR") << "Cannot use chunking in window mode!" << std::endl;
			return false;
		}

		std::cerr << utility::timestamp("LOG","READER") << "Opening " << settings.in << "..." << std::endl;

		twk_reader reader;
		if(reader.Open(settings.in) == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to open file: " << settings.in << "..." << std::endl;
			return false;
		}

		std::cerr << utility::timestamp("LOG") << "Samples: " << utility::ToPrettyString(reader.hdr.GetNumberSamples()) << "..." << std::endl;

		twk1_blk_iterator bit;
		bit.stream = reader.stream;

		Timer timer;

		//settings.ldd_load_type = TWK_LDD_ALL;
		if(settings.low_memory && settings.force_phased && settings.bitmaps){
			settings.ldd_load_type = TWK_LDD_BITMAP;
		} else if(settings.low_memory && settings.force_phased){
			settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
		} else if(settings.force_phased){
			settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
		} else if(settings.forced_unphased){
			settings.ldd_load_type = TWK_LDD_VEC;
		} else {
			settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
		}

		if(settings.ival_strings.size() == 0) n_blks = reader.index.n;
		else {
			if(this->BuildIntervals(reader, bit) == false)
				return false;
		}

		if(settings.window) settings.c_chunk = 0;

		twk_ld_balancer balancer;
		if(balancer.Build(this->n_blks, settings.n_chunks, settings.c_chunk) == false){
			return false;
		}

		std::cerr << utility::timestamp("LOG","BALANCING") << "Using ranges [" << balancer.fromL << "-" << balancer.toL << "," << balancer.fromR << "-" << balancer.toR << "] in " << (settings.window ? "window mode" : "square mode") <<"..." << std::endl;

		if(this->n_blks == 0){
			std::cerr << utility::timestamp("ERROR") << "No valid data available..." << std::endl;
			return true;
		}

		uint32_t n_variants = 0;
		if(balancer.diag){
			for(int i = balancer.fromL; i < balancer.toL; ++i) n_variants += reader.index.ent[i].n;
		} else {
			for(int i = balancer.fromL; i < balancer.toL; ++i) n_variants += reader.index.ent[i].n;
			for(int i = balancer.fromR; i < balancer.toR; ++i) n_variants += reader.index.ent[i].n;
		}

		std::cerr << utility::timestamp("LOG") << utility::ToPrettyString(n_variants) << " variants from " << utility::ToPrettyString(balancer.n_m) << " blocks..." << std::endl;
		std::cerr << utility::timestamp("LOG","PARAMS") << settings.GetString() << std::endl;
		std::cerr << utility::timestamp("LOG") << "Performing: " << utility::ToPrettyString(((uint64_t)n_variants * n_variants - n_variants) / 2) << " variant comparisons..." << std::endl;

		twk_ld_progress progress;
		progress.n_s = reader.hdr.GetNumberSamples();
		if(settings.window == false){
			if(balancer.diag)
				progress.n_cmps = ((uint64_t)n_variants * n_variants - n_variants) / 2;
			else {

			}
		}

		// Append literal string.
		std::string calc_string = "\n##tomahawk_calcVersion=" + std::to_string(VERSION) + "\n";
		calc_string += "##tomahawk_calcCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + utility::datetime() + "\n";
		reader.hdr.literals_ += calc_string;

		// New index
		IndexOutput index(reader.hdr.GetNumberContigs());

		twk_ld_engine::ep f; // function pointer array
		f[0] = &twk_ld_engine::PhasedVectorized;
		f[1] = &twk_ld_engine::PhasedVectorizedNoMissing;
		f[2] = &twk_ld_engine::UnphasedVectorized;
		f[3] = &twk_ld_engine::UnphasedVectorizedNoMissing;
		f[4] = &twk_ld_engine::PhasedRunlength;
		f[5] = &twk_ld_engine::PhasedList;
		f[6] = &twk_ld_engine::UnphasedRunlength;
		f[7] = &twk_ld_engine::PhasedBitmap;
		f[8] = &twk_ld_engine::PhasedListSpecial;
		f[9] = &twk_ld_engine::PhasedListVector;

		const std::vector<std::string> method_names = {"PhasedVectorized",
				"PhasedVectorizedNoMissing",
				"UnphasedVectorized",
				"UnphasedVectorizedNoMissing",
				"PhasedRunlength",
				"PhasedList",
				"UnphasedRunlength",
				"PhasedBitmap",
				"PhasedListSpecial",
				"PhasedListBV"};

		uint8_t unpack_list[10];
		memset(unpack_list, TWK_LDD_VEC, 10);
		unpack_list[4] = TWK_LDD_NONE;
		unpack_list[5] = TWK_LDD_LIST;
		unpack_list[6] = TWK_LDD_NONE;
		unpack_list[7] = TWK_LDD_BITMAP;
		unpack_list[8] = TWK_LDD_VEC | TWK_LDD_LIST;
		unpack_list[9] = TWK_LDD_VEC | TWK_LDD_LIST;

		twk_ld_perf perfs[10];
		uint32_t perf_size = reader.hdr.GetNumberSamples()*4 + 2;
		for(int i = 0; i < 10; ++i){
			perfs[i].cycles = new uint64_t[perf_size];
			perfs[i].freq   = new uint64_t[perf_size];
			memset(perfs[i].cycles, 0, perf_size*sizeof(uint64_t));
			memset(perfs[i].freq,   0, perf_size*sizeof(uint64_t));
		}

		//prepare_shuffling_dictionary();
		for(int method = 8; method < 10; ++method){
			if(method == 2 || method == 3 || method == 4 || method == 6) continue;
			std::cerr << utility::timestamp("LOG","PROGRESS") << "Starting method: " << method_names[method] << "..." << std::endl;
			settings.ldd_load_type = unpack_list[method];
			if(this->LoadBlocks(reader, bit, balancer, settings) == false){
				return false;
			}

			twk_ld_ticker ticker;
			ticker = balancer;
			ticker = settings;
			ticker.ldd  = ldd;

			timer.Start();
			twk_ld_slave s;
			s.ldd    = ldd;
			s.n_s    = reader.hdr.GetNumberSamples();
			s.ticker = &ticker;
			s.engine.SetSamples(reader.hdr.GetNumberSamples());
			s.engine.SetBlocksize(settings.b_size);
			s.engine.progress = &progress;
			//s.engine.writer   = writer;
			s.engine.index    = &index;
			s.engine.settings = settings;
			s.progress = &progress;
			s.settings = &settings;

			s.CalculatePerformance(f[method], &perfs[method]);

			std::cerr << utility::timestamp("LOG","PROGRESS") << "Finished method: " << method_names[method] << ". Elapsed time=" << timer.ElapsedString() << std::endl;
			//std::cerr << "m=" << method << " " << utility::ToPrettyString(1) << " " << utility::ToPrettyString((uint64_t)((float)1/timer.Elapsed().count())) << "/s " << utility::ToPrettyString((uint64_t)((float)1*reader.hdr.GetNumberSamples()/timer.Elapsed().count())) << "/gt/s " << timer.ElapsedString() << " " << utility::ToPrettyString(1) << std::endl;
		}

		for(int i = 0; i < 10; ++i){
			for(int j = 0; j < perf_size; ++j){
				// Print if non-zero
				if(perfs[i].freq[j] != 0)
					std::cout << i << "\t" << j << "\t" << (double)perfs[i].cycles[j]/perfs[i].freq[j] << "\t" << perfs[i].freq[j] << "\t" << perfs[i].cycles[j] << '\n';
			}
			delete[] perfs[i].cycles;
		}

		return true;
	}

public:
	uint32_t n_blks, m_blks, n_vnts, n_tree;
	twk1_ldd_blk* ldd;
	twk1_block_t* ldd2;
	twk_ld_settings settings;
	twk_intervals_twk intervals;
};

}

#endif
