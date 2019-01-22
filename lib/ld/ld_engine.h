#ifndef TWK_LD_ENGINE_H_
#define TWK_LD_ENGINE_H_

#include "core.h"
#include "twk_reader.h"
#include "writer.h"

// Make sure they are not in the API
#include "ld/ld_structs.h"
#include "ld/ld_progress.h"
#include "ld/ld_balancing.h"
#include "ld/ld_unpacker.h"

namespace tomahawk {

// Method 0: None
// Method 1: Performance + no math
// Method 2: Debug correctness
// Method 3: Print LD difference Phased and Unphased
#define TWK_SLAVE_DEBUG_MODE 0

// Accessors for ease-of-use
#define TWK_LD_SIMD_REFREF 0
#define TWK_LD_SIMD_ALTREF 1
#define TWK_LD_SIMD_REFALT 2
#define TWK_LD_SIMD_ALTALT 3
#define TWK_LD_REFREF 0
#define TWK_LD_ALTREF 1
#define TWK_LD_REFALT 4
#define TWK_LD_ALTALT 5

// Parameter thresholds for FLAGs
#define TWK_LOW_AC_THRESHOLD        5
#define TWK_INVALID_HWE_THRESHOLD   1e-4
#define TWK_LONG_RANGE_THRESHOLD    500e3
#define TWK_MINIMUM_ALLOWED_ALLELES 5		// Minimum number of alleles required for math to work out in the unphased case
#define TWK_ALLOWED_ROUNDING_ERROR  0.00001

#define TWK_HAP_FREQ(A,POS) ((double)(A).alleleCounts[POS] / (A).totalHaplotypeCounts)

#if SLAVE_DEBUG_MODE == 3
static uint32_t twk_debug_pos1   = 0;
static uint32_t twk_debug_pos1_2 = 0;
static uint32_t twk_debug_pos2   = 0;
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

#define UNPHASED_UPPER_MASK     (uint64_t)170  // 10101010b
#define UNPHASED_LOWER_MASK     (uint64_t)85   // 01010101b
#define UNPHASED_UPPER_MASK_64  ((UNPHASED_UPPER_MASK << 56) | (UNPHASED_UPPER_MASK << 48) | (UNPHASED_UPPER_MASK << 40) | (UNPHASED_UPPER_MASK << 32) | (UNPHASED_UPPER_MASK << 24) | (UNPHASED_UPPER_MASK << 16) | (UNPHASED_UPPER_MASK << 8) | (UNPHASED_UPPER_MASK))  // 10101010 10101010 10101010 10101010b
#define UNPHASED_LOWER_MASK_64  ((UNPHASED_LOWER_MASK << 56) | (UNPHASED_LOWER_MASK << 48) | (UNPHASED_LOWER_MASK << 40) | (UNPHASED_LOWER_MASK << 32) | (UNPHASED_LOWER_MASK << 24) | (UNPHASED_LOWER_MASK << 16) | (UNPHASED_LOWER_MASK << 8) | (UNPHASED_LOWER_MASK))  // 01010101 01010101 01010101 01010101b
#define FILTER_UNPHASED_64(A, B)            ((((((A) & UNPHASED_UPPER_MASK_64) | ((B) & UNPHASED_LOWER_MASK_64)) & UNPHASED_LOWER_MASK_64) << 1) & (A))
#define FILTER_UNPHASED_64_PAIR(A, B, C, D) ((FILTER_UNPHASED_64((A), (B)) >> 1) | FILTER_UNPHASED_64((C), (D)))
#define FILTER_UNPHASED_64_SPECIAL(A)       ((((A) >> 1) & (A)) & UNPHASED_LOWER_MASK_64)

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
 * Supportive counter struct for when using the vectorized instructions
 * for computing linkage-disequilibrium.
 */
struct twk_ld_simd {
public:
	twk_ld_simd(void);
	~twk_ld_simd();

public:
	uint64_t *counters;
} __attribute__((aligned(16)));

/**<
 * Generic supportive counter struct for when calculating linkage-disequilibrium.
 * This struct have internal counters for haplotype and allele counts and for the
 * total number of allele/haplotype counts.
 */
struct twk_ld_count {
	twk_ld_count();
	~twk_ld_count();

	void ResetPhased(void);
	void ResetUnphased(void);

	// Counters
	uint64_t alleleCounts[171]; // equivalent to 1010,1010b = (MISS,MISS),(MISS,MISS)
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
	twk_ld_engine();
	~twk_ld_engine();
	twk_ld_engine(const twk_ld_engine& other) = delete;
	twk_ld_engine& operator=(const twk_ld_engine& other) = delete;
	twk_ld_engine(twk_ld_engine&& other) = delete;
	twk_ld_engine& operator=(twk_ld_engine&& other) = delete;

	/**<
	 * Set the number of samples in the target file. This function is mandatory
	 * to call prior to computing as it sets essential paramters such as the number
	 * of samples, the byte alignment of vectorized bit-vectors and the adjustment
	 * given the unused overhangs.
	 * @param samples Number of samples in the target file.
	 */
	void SetSamples(const uint32_t samples);
	void SetBlocksize(const uint32_t s);

	/**<
	 * Phased/unphased functions for calculating linkage-disequilibrium. These
	 * functions are all prefixed with Phased_ or Unphased_. All these functions
	 * share the same interface in order to allow them being targetted by a shared
	 * function pointer.
	 *
	 * Algorithms:
	 *    *RunLength  compared two word-aligned RLE objects in O(|A| + |B| + 1) time.
	 *    *Vectorized use machine-optimized SIMD instructions to horizontally compare
	 *                register-aligned arrays.
	 *    *List       combines an inverted index mapping alt-positions to a word
	 *                offset within the bitvector used in *Vectorized functions. These
	 *                functions use select operations on the bitvector directly.
	 *    *ListVector This variation maps to register offsets directly followed by a
	 *                popcount operation instead of individual select queries. This
	 *                effectively converge into *Vectorize functions when the size of
	 *                the inverted index equals that of the number of registers.
	 *
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
	//bool PhasedListSpecial(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool PhasedVectorized(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool PhasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool PhasedBitmap(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool PhasedMath(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2);
	bool UnphasedRunlength(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool UnphasedVectorized(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool UnphasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);
	bool UnphasedList(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf = nullptr);

	// Unphased math.
	bool UnphasedMath(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2);
	double ChiSquaredUnphasedTable(const double target, const double p, const double q);
	bool ChooseF11Calculate(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2,const double target, const double p, const double q);

	/**<
	 * Compress the internal output buffers consisting of twk1_two_t records.
	 * There are two separate buffers, Fwd and Rev, corresponding to the upper
	 * and lower triangular of output associations. This way output data is
	 * more partially sorted compared to interleaving these two streams.
	 *
	 * The CompressBlock function encapsulates both Fwd and Rev compress
	 * functions. Calling this function is sufficient and recommended.
	 *
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
	uint32_t phased_unbalanced_adjustment; // Remainder in math
	uint32_t unphased_unbalanced_adjustment; // Remainder in math
	uint64_t t_out; // number of bytes written
	uint64_t n_method[10]; // Number of times a tgt function was used

	uint64_t* mask_placeholder; // placeholder all-0 mask.

	IndexEntryOutput irecF, irecR;
	ZSTDCodec zcodec; // reusable zstd codec instance with internal context.

	twk_ld_settings settings;
	twk1_two_block_t blk_f, blk_r;
	twk_buffer_t ibuf, obuf; // buffers
	twk_ld_simd  helper_simd; // simd counters
	twk_ld_count helper; // regular counters
	twk1_two_t cur_rcd;
	IndexOutput* index;
	twk_writer_t* writer;
	twk_ld_progress* progress;
	uint32_t* list_out;
};

/**<
 * Thread slave for computing linkage-disequilibrium in parallel.
 */
struct twk_ld_slave {
	twk_ld_slave();
	~twk_ld_slave();

	/**<
	 * Primary subroutine for starting the slave thread to compute its designated
	 * region(s) of linkage-disequilibrium.
	 * @return Returns a pointer to the spawned slave thread.
	 */
	std::thread* Start();

	/**<
	 * Update local blocks. Either copies pointer to pre-allocated containers or
	 * construct new on-the-fly given the memory parametrisation. Internally
	 * ivokes the subroutines prefixed with `UpdateBlocks`.
	 * @param blks  Dst twk1_ldd_blk array of size 2.
	 * @param from  Integer start offset.
	 * @param to    Interger end offset.
	 */
	void UpdateBlocks(twk1_ldd_blk* blks, const uint32_t& from, const uint32_t& to);

	/**<
	 * Subroutine for retrieving a pair of twk1_ldd_blk from the pre-loaded arrays.
	 * @param blks Dst pointers to twk1_ldd_blks.
	 * @param from Virtual offset for A.
	 * @param to   Virtual offset for B.
	 */
	void UpdateBlocksPreloaded(twk1_ldd_blk* blks, const uint32_t& from, const uint32_t& to);

	/**<
	 * Subroutine for online pre-computing of a pair of twk1_ldd_blk. This
	 * subroutine is used exclusively in low-memory mode.
	 * @param blks Dst pointers to twk1_ldd_blks.
	 * @param from Virtual offset for A.
	 * @param to   Virtual offset for B.
	 */
	void UpdateBlocksGenerate(twk1_ldd_blk* blks, const uint32_t& from, const uint32_t& to);

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
	            twk_ld_perf* perf = nullptr);

	void Unphased(const twk1_t* rcds0,
	              const twk1_t* rcds1,
	              const twk1_ldd_blk* blocks,
	              const uint8_t type,
	              twk_ld_perf* perf = nullptr);

	bool CalculatePhased(twk_ld_perf* perf = nullptr);
	bool CalculateUnphased(twk_ld_perf* perf = nullptr);
	bool CalculatePhasedBitmap(twk_ld_perf* perf = nullptr);
	bool CalculatePhasedBitmapWindow(twk_ld_perf* perf = nullptr);
	bool CalculatePhasedWindow(twk_ld_perf* perf = nullptr);
	bool CalculateUnphasedWindow(twk_ld_perf* perf = nullptr);
	bool Calculate(twk_ld_perf* perf = nullptr);
	bool CalculatePerformance(twk_ld_engine::func f, twk_ld_perf* perf = nullptr);
	bool CalculateSingle(twk_ld_perf* perf = nullptr);

public:
	uint32_t n_s, n_total;
	uint32_t i_start, j_start, prev_i, prev_j, n_cycles;

	twk_ld_dynamic_balancer* ticker;
	std::thread* thread;
	twk1_ldd_blk* ldd;
	twk_ld_progress* progress;
	twk_ld_settings* settings;
	twk_ld_engine engine;
};

}


#endif /* LIB_LD_LD_ENGINE_H_ */
