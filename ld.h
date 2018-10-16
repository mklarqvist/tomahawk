#ifndef TWK_LD_H_
#define TWK_LD_H_

#include <cassert>
#include "core.h"

// Method 1: None: Input-specified (default)
// Method 2: Phased Vectorized No-Missing
// Method 3: Count comparisons for A1
// Method 4: Unphased A1 and unphased A2
// Method 5: Phased A1 and Phased A2
// Method 6: All algorithms comparison (debug)
// Method 7: All algorithms run-time output (debug)
#define SLAVE_DEBUG_MODE 1

namespace tomahawk{

// Parameter thresholds for FLAGs
#define LOW_MAF_THRESHOLD       0.01
#define LOW_HWE_THRESHOLD       1e-6
#define LONG_RANGE_THRESHOLD    500e3
#define MINIMUM_ALLOWED_ALLELES 5		// Minimum number of alleles required for math to work out in the unphased case

// SIMD trigger
#if SIMD_AVAILABLE == 1

#ifdef _popcnt64
#define POPCOUNT_ITER	_popcnt64
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

#define POPCOUNT(A, B) {									\
	__m256i tempA = _mm512_extracti64x4_epi64(B, 0);		\
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 0));		\
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 1)); 	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 2)); 	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 3)); 	\
	tempA = _mm512_extracti64x4_epi64(B, 1); 				\
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 0));		\
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 1)); 	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 2)); 	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(tempA, 3)); 	\
}

#define FILTER_UNPHASED(A, B)			 _mm512_and_si512(_mm512_slli_epi64(_mm512_and_si512(_mm512_or_si512(_mm512_and_si512(A, maskUnphasedHigh),_mm512_and_si512(B, maskUnphasedLow)), maskUnphasedLow), 1), A)
#define FILTER_UNPHASED_PAIR(A, B, C, D) _mm512_or_si512(_mm512_srli_epi64(FILTER_UNPHASED(A, B), 1), FILTER_UNPHASED(C, D))
#define FILTER_UNPHASED_SPECIAL(A)		 _mm512_and_si512(_mm512_and_si512(_mm512_srli_epi64(A, 1), A), maskUnphasedLow)


#elif SIMD_VERSION == 5 // AVX2
#define VECTOR_TYPE	__m256i
const VECTOR_TYPE ONE_MASK         = _mm256_set1_epi8(255); // 11111111b
const VECTOR_TYPE maskUnphasedHigh = _mm256_set1_epi8(UNPHASED_UPPER_MASK);	// 10101010b
const VECTOR_TYPE maskUnphasedLow  = _mm256_set1_epi8(UNPHASED_LOWER_MASK);	// 01010101b

#define PHASED_ALTALT(A,B)	        _mm256_and_si256(A, B)
#define PHASED_REFREF(A,B)	        _mm256_and_si256(_mm256_xor_si256(A, ONE_MASK), _mm256_xor_si256(B, ONE_MASK))
#define PHASED_ALTREF(A,B)	        _mm256_and_si256(_mm256_xor_si256(A, B), B)
#define PHASED_REFALT(A,B)	        _mm256_and_si256(_mm256_xor_si256(A, B), A)
#define PHASED_ALTALT_MASK(A,B,M)	_mm256_and_si256(PHASED_ALTALT(A, B), M)
#define PHASED_REFREF_MASK(A,B,M)	_mm256_and_si256(PHASED_REFREF(A, B), M)
#define PHASED_ALTREF_MASK(A,B,M)	_mm256_and_si256(PHASED_ALTREF(A, B), M)
#define PHASED_REFALT_MASK(A,B,M)	_mm256_and_si256(PHASED_REFALT(A, B), M)
#define MASK_MERGE(A,B)		        _mm256_xor_si256(_mm256_or_si256(A, B), ONE_MASK)

// Software intrinsic popcount
#define POPCOUNT(A, B) {							\
	A += POPCOUNT_ITER(_mm256_extract_epi64(B, 0));	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(B, 1));	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(B, 2));	\
	A += POPCOUNT_ITER(_mm256_extract_epi64(B, 3));	\
}

#define FILTER_UNPHASED(A, B)			 _mm256_and_si256(_mm256_slli_epi64(_mm256_and_si256(_mm256_or_si256(_mm256_and_si256(A, maskUnphasedHigh),_mm256_and_si256(B, maskUnphasedLow)), maskUnphasedLow), 1), A)
#define FILTER_UNPHASED_PAIR(A, B, C, D) _mm256_or_si256(_mm256_srli_epi64(FILTER_UNPHASED(A, B), 1), FILTER_UNPHASED(C, D))
#define FILTER_UNPHASED_SPECIAL(A)		 _mm256_and_si256(_mm256_and_si256(_mm256_srli_epi64(A, 1), A), maskUnphasedLow)

#elif SIMD_VERSION >= 2 // SSE2+
#define VECTOR_TYPE	__m128i
const VECTOR_TYPE ONE_MASK         = _mm_set1_epi8(255); // 11111111b
const VECTOR_TYPE maskUnphasedHigh = _mm_set1_epi8(UNPHASED_UPPER_MASK);	// 10101010b
const VECTOR_TYPE maskUnphasedLow  = _mm_set1_epi8(UNPHASED_LOWER_MASK);	// 01010101b

#define PHASED_ALTALT(A,B)	        _mm_and_si128(A, B)
#define PHASED_REFREF(A,B)	        _mm_and_si128(_mm_xor_si128(A, ONE_MASK), _mm_xor_si128(B, ONE_MASK))
#define PHASED_ALTREF(A,B)	        _mm_and_si128(_mm_xor_si128(A, B), B)
#define PHASED_REFALT(A,B)	        _mm_and_si128(_mm_xor_si128(A, B), A)
#define PHASED_ALTALT_MASK(A,B,M)	_mm_and_si128(PHASED_ALTALT(A, B), M)
#define PHASED_REFREF_MASK(A,B,M)	_mm_and_si128(PHASED_REFREF(A, B), M)
#define PHASED_ALTREF_MASK(A,B,M)	_mm_and_si128(PHASED_ALTREF(A, B), M)
#define PHASED_REFALT_MASK(A,B,M)	_mm_and_si128(PHASED_REFALT(A, B), M)
#define MASK_MERGE(A,B)		        _mm_xor_si128(_mm_or_si128(A, B), ONE_MASK)

#if SIMD_VERSION >= 3
#define POPCOUNT(A, B) {								\
	A += POPCOUNT_ITER(_mm_extract_epi64(B, 0));		\
	A += POPCOUNT_ITER(_mm_extract_epi64(B, 1));		\
}
#else
#define POPCOUNT(A, B) { \
	uint64_t temp = _mm_extract_epi16(B, 0) << 6 | _mm_extract_epi16(B, 1) << 4 | _mm_extract_epi16(B, 2) << 2 | _mm_extract_epi16(B, 3); \
	A += POPCOUNT_ITER(temp);		\
	temp = _mm_extract_epi16(B, 4) << 6 | _mm_extract_epi16(B, 5) << 4 | _mm_extract_epi16(B, 6) << 2 | _mm_extract_epi16(B, 7); \
	A += POPCOUNT_ITER(temp);		\
}
#endif

#define FILTER_UNPHASED(A, B)			 _mm_and_si128(_mm_slli_epi64(_mm_and_si128(_mm_or_si128(_mm_and_si128(A, maskUnphasedHigh),_mm_and_si128(B, maskUnphasedLow)), maskUnphasedLow), 1), A)
#define FILTER_UNPHASED_PAIR(A, B, C, D) _mm_or_si128(_mm_srli_epi64(FILTER_UNPHASED(A, B), 1), FILTER_UNPHASED(C, D))
#define FILTER_UNPHASED_SPECIAL(A)		 _mm_and_si128(_mm_and_si128(_mm_srli_epi64(A, 1), A), maskUnphasedLow)
#endif
#endif // ENDIF SIMD_AVAILABLE == 1

struct twk_ld_simd {
public:
	twk_ld_simd(void);
	~twk_ld_simd();

public:
	uint64_t *counters;
	uint8_t  *scalarA, *scalarB, *scalarC, *scalarD;
} __attribute__((aligned(16)));

struct twk_ld_count {
	twk_ld_count(){}
	~twk_ld_count(){}

	void resetPhased(void){
		this->alleleCounts[0] = 0;
		this->alleleCounts[1] = 0;
		this->alleleCounts[4] = 0;
		this->alleleCounts[5] = 0;
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
	double alleleCounts[171];
	double haplotypeCounts[4];
};

class LDEngine {
public:

void SetSamples(const uint32_t samples){
	n_samples  = samples;
	byte_width = ceil((double)samples/4);
	byte_aligned_end = byte_width/(GENOTYPE_TRIP_COUNT/4)*(GENOTYPE_TRIP_COUNT/4);
	vector_cycles    = byte_aligned_end*4/GENOTYPE_TRIP_COUNT;
	phased_unbalanced_adjustment   = (samples*2)%8;
	unphased_unbalanced_adjustment = samples%4;
}

bool CalculateLDPhasedSimple(const twk_igt_list& ref, const twk_igt_list& tgt){
	helper.resetPhased();

	const uint32_t n_cycles = ref.l_list < tgt.l_list ? ref.l_list : tgt.l_list;
	const uint32_t n_total  = ref.l_list + tgt.l_list;
	uint32_t n_same = 0;

	if(ref.l_list >= tgt.l_list){
		for(uint32_t i = 0; i < n_cycles; ++i) n_same += ref.get(tgt.list[i]);
	} else {
		for(uint32_t i = 0; i < n_cycles; ++i) n_same += tgt.get(ref.list[i]);
	}

	helper.haplotypeCounts[0] = 2*n_samples - (n_total - n_same);
	//std::cerr << 2*this->n_samples - (n_total - n_same) << '\n';
	//assert(2*n_samples - (n_total - n_same) < 2*n_samples);
	//std::cerr << helper.haplotypeCounts[0] << " ";

	//return(this->CalculateLDPhasedMathSimple(block1, block2));
	return true;
}

bool PhasedVectorizedNoMissingNoTable(const twk_igt_vec& block1, const twk_igt_vec& block2){
	helper.resetPhased();
	helper_simd.counters[0] = 0;

	const uint8_t* const arrayA = block1.data;
	const uint8_t* const arrayB = block2.data;

#if SIMD_AVAILABLE == 1
	const uint32_t frontSmallest = block1.front_zero < block2.front_zero ? block1.front_zero : block2.front_zero;
	const uint32_t tailSmallest  = block1.tail_zero  < block2.tail_zero  ? block1.tail_zero  : block2.tail_zero;
	const uint32_t frontBonus    = block1.front_zero != frontSmallest ? block1.front_zero : block2.front_zero;
	const uint32_t tailBonus     = block1.tail_zero  != tailSmallest  ? block1.tail_zero  : block2.tail_zero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	VECTOR_TYPE __intermediate;

#define ITER_SHORT {											\
	__intermediate  = PHASED_ALTALT(vectorA[i], vectorB[i]);	\
	POPCOUNT(helper_simd.counters[0], __intermediate);			\
	i += 1;														\
}

	uint32_t i = frontSmallest;
	for( ; i < frontBonus; )                         ITER_SHORT
	for( ; i < this->vector_cycles - tailBonus; )    ITER_SHORT
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT

#undef ITER_SHORT
	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif

	//std::cerr << "k=" << k << "/" << this->byte_width << std::endl;
	for(; k+8 < this->byte_width; k += 8){
		for(uint32_t l = 0; l < 8; ++l)
			helper_simd.scalarB[l] = arrayA[k+l] & arrayB[k+l];
		helper_simd.counters[0] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarA));
	}

	for(; k < this->byte_width; ++k){
		helper_simd.counters[0] += POPCOUNT_ITER(((arrayA[k] & arrayB[k]) & 255));
	}
	//helper.haplotypeCounts[0] = (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 + helper_simd.counters[0] - this->phased_unbalanced_adjustment;
	//std::cerr << helper.haplotypeCounts[0] << std::endl;

	//this->setFLAGs(block1, block2);
	//return(this->CalculateLDPhasedMathSimple(block1, block2));
	//exit(1);
	return true;
}

/*
bool MathPhasedSimple(const block_type& block1, const block_type& block2){
	//++this->possible;

	// If haplotype count for (0,0) >
	//if(helper[0] > 2*n_samples - this->parameters.minimum_sum_alternative_haplotype_count)
	//	return false;

	// D = (joint HOM_HOM) - (HOM_A * HOM_B) = pAB - pApB
	const double divisor = block1.currentMeta().AF * (1 - block1.currentMeta().AF) * block2.currentMeta().AF * (1 - block2.currentMeta().AF);
	if(divisor == 0) return false;

	double D, Dprime, R2, R;

	D  = helper.haplotypeCounts[0]/(2*n_samples) - block1.currentMeta().AF*block2.currentMeta().AF;
	R2 = D*D / divisor ;
	R  = sqrt(R2);

	if(D < 0){
		if( block1.currentMeta().AF * block2.currentMeta().AF < (1 - block1.currentMeta().AF) * (1 - block2.currentMeta().AF)){
			Dprime = D / (block1.currentMeta().AF * block2.currentMeta().AF);
		} else {
			Dprime = D / ( (1 - block1.currentMeta().AF) * (1 - block2.currentMeta().AF) );
		}
	} else { // D >= 0
		if((1-block1.currentMeta().AF) * block2.currentMeta().AF < block1.currentMeta().AF * (1 - block2.currentMeta().AF)){
			Dprime = D / ((1 - block1.currentMeta().AF) * block2.currentMeta().AF);
		} else {
			Dprime = D / (block1.currentMeta().AF * (1 - block2.currentMeta().AF));
		}
	}

	return true;

}
*/

public:
	uint32_t n_samples;
	uint32_t byte_width; // Number of bytes required per variant site
	uint32_t byte_aligned_end; // End byte position
	uint32_t vector_cycles; // Number of SIMD cycles (genotypes/2/vector width)
	uint32_t phased_unbalanced_adjustment; // Modulus remainder
	uint32_t unphased_unbalanced_adjustment; // Modulus remainder

	twk_ld_simd  helper_simd;
	twk_ld_count helper;
};

}


#endif /* LD_H_ */
