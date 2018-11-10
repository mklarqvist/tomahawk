#include <chrono>
#include <cstdlib>

#include "ld.h"
#include "fisher_math.h"

namespace tomahawk {

twk_ld_simd::twk_ld_simd(void) :
#if SIMD_AVAILABLE == 1
	counters((uint64_t*)_mm_malloc(sizeof(uint64_t)*16, 16)),
	scalarA((uint8_t*)_mm_malloc(sizeof(uint8_t)*8, 16)),
	scalarB((uint8_t*)_mm_malloc(sizeof(uint8_t)*8, 16)),
	scalarC((uint8_t*)_mm_malloc(sizeof(uint8_t)*8, 16)),
	scalarD((uint8_t*)_mm_malloc(sizeof(uint8_t)*8, 16))
#else
	counters(new uint64_t[16]),
	scalarA(new uint8_t[8]),
	scalarB(new uint8_t[8]),
	scalarC(new uint8_t[8]),
	scalarD(new uint8_t[8])
#endif
{
	memset(this->counters, 0, sizeof(uint64_t)*16);
}

twk_ld_simd::~twk_ld_simd(){
#if SIMD_AVAILABLE == 1
	_mm_free(this->counters);
	_mm_free(this->scalarA);
	_mm_free(this->scalarB);
	_mm_free(this->scalarC);
	_mm_free(this->scalarD);
#else
	delete [] this->counters;
	delete [] this->scalarA;
	delete [] this->scalarB;
	delete [] this->scalarC;
	delete [] this->scalarD;
#endif
}


//

bool twk_ld_engine::PhasedList(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, twk_ld_perf* perf){
	helper.resetPhased();
	const twk_igt_list& ref = b1.list[p1];
	const twk_igt_list& tgt = b2.list[p2];
	const uint32_t n_cycles = ref.l_list < tgt.l_list ? ref.l_list : tgt.l_list;
	//const uint32_t n_total  = ref.l_list + tgt.l_list;
	//uint32_t n_same = 0;

	// Debug timings
	#if SLAVE_DEBUG_MODE == 1
		typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
		auto t0 = std::chrono::high_resolution_clock::now();
	#endif

	if(ref.l_list >= tgt.l_list){
		for(uint32_t i = 0; i < n_cycles; ++i) helper.alleleCounts[5] += ref.get(tgt.list[i]);
	} else {
		for(uint32_t i = 0; i < n_cycles; ++i) helper.alleleCounts[5] += tgt.get(ref.list[i]);
	}

	helper.alleleCounts[4] = ref.l_list - helper.alleleCounts[5];
	helper.alleleCounts[1] = tgt.l_list - helper.alleleCounts[5];
	helper.alleleCounts[0] = 2*n_samples - ((ref.l_list + tgt.l_list) - helper.alleleCounts[5]);
	//std::cerr << "phased-list=" << ref.l_list << " and " << tgt.l_list << std::endl;
	//std::cerr << 2*this->n_samples - (n_total - n_same) << '\n';
	//assert(2*n_samples - (n_total - n_same) < 2*n_samples);
	//std::cerr << helper.haplotypeCounts[0] << " ";

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[n_cycles] += ticks_per_iter.count();
	++perf->freq[n_cycles];
#endif

#if SLAVE_DEBUG_MODE == 2
	std::cerr << "list=" << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[4] << "," << helper.alleleCounts[5] << std::endl;
#endif

#if SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::PhasedBitmap(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, twk_ld_perf* perf){
	helper.resetPhased();
	const twk1_ldd_blk::bitmap_type& refB = b1.bitmap[p1];
	const twk1_ldd_blk::bitmap_type& tgtB = b2.bitmap[p2];
	helper.alleleCounts[5] = refB.logicalandcount(tgtB);
	helper.alleleCounts[4] = b1.blk->rcds[p1].ac - helper.alleleCounts[5];
	helper.alleleCounts[1] = b2.blk->rcds[p2].ac - helper.alleleCounts[5];
	helper.alleleCounts[0] = 2*n_samples - ((b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac) - helper.alleleCounts[5]);

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[n_cycles] += ticks_per_iter.count();
	++perf->freq[n_cycles];
#endif

#if SLAVE_DEBUG_MODE == 2
	std::cerr << "bitmap=" << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[4] << "," << helper.alleleCounts[5] << std::endl;
#endif

#if SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::PhasedVectorized(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, twk_ld_perf* perf){
#if SLAVE_DEBUG_MODE == 0
	if(b1.blk->rcds[p1].gt_missing == false && b2.blk->rcds[p2].gt_missing == false){
		return(this->PhasedVectorizedNoMissing(b1,p1,b2,p2,perf));
		//return(this->CalculateLDPhasedVectorizedNoMissingNoTable(helper, block1, block2));
	}
#endif

	helper.resetPhased();
	helper_simd.counters[0] = 0;
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;
	// Data
	const twk_igt_vec& block1 = b1.vec[p1];
	const twk_igt_vec& block2 = b2.vec[p2];
	const uint8_t* const arrayA = (const uint8_t* const)block1.data;
	const uint8_t* const arrayB = (const uint8_t* const)block2.data;
	const uint8_t* const arrayA_mask = (const uint8_t* const)block1.mask;
	const uint8_t* const arrayB_mask = (const uint8_t* const)block2.mask;

#if SIMD_AVAILABLE == 1
	const uint32_t frontSmallest = block1.front_zero < block2.front_zero ? block1.front_zero : block2.front_zero;
	const uint32_t tailSmallest  = block1.tail_zero  < block2.tail_zero  ? block1.tail_zero  : block2.tail_zero;
	const uint32_t frontBonus    = block1.front_zero != frontSmallest ? block1.front_zero : block2.front_zero;
	const uint32_t tailBonus     = block1.tail_zero  != tailSmallest  ? block1.tail_zero  : block2.tail_zero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	const VECTOR_TYPE* const vectorA_mask = (const VECTOR_TYPE* const)arrayA_mask;
	const VECTOR_TYPE* const vectorB_mask = (const VECTOR_TYPE* const)arrayB_mask;

	VECTOR_TYPE __intermediate, masks;

	uint32_t i = frontSmallest;

// Debug timings
#if SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_SHORT {                                                     \
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);              \
	__intermediate  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_REFREF], __intermediate);       \
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_ALTREF], __intermediate);       \
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_REFALT], __intermediate);       \
	i += 1;                                                              \
}

#define ITER {                                                           \
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);              \
	__intermediate  = PHASED_ALTALT_MASK(vectorA[i], vectorB[i], masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_ALTALT], __intermediate);       \
	__intermediate  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_REFREF], __intermediate);       \
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_ALTREF], __intermediate);       \
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_REFALT], __intermediate);       \
	i += 1;                                                              \
}

	for( ; i < frontBonus; ) 					  	 ITER_SHORT // Not possible to be ALT-ALT
	for( ; i < this->vector_cycles - tailBonus; )  	 ITER
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT // Not possible to be ALT-ALT

#undef ITER
#undef ITER_SHORT
	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif

	uint8_t mask;
	for(; k+8 < this->byte_width; k += 8){
		for(uint32_t l = 0; l < 8; ++l){
			mask = ~(arrayA_mask[k+l] | arrayB_mask[k+l]);
			helper_simd.scalarA[l] = (arrayA[k+l] & arrayB[k+l]) & mask;
			helper_simd.scalarB[l] = ((~arrayA[k+l]) & (~arrayB[k+l])) & mask;
			helper_simd.scalarC[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayA[k+l]) & mask;
			helper_simd.scalarD[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayB[k+l]) & mask;
		}
		helper_simd.counters[TWK_LD_SIMD_REFREF] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarB));
		helper_simd.counters[TWK_LD_SIMD_REFALT] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarC));
		helper_simd.counters[TWK_LD_SIMD_ALTREF] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarD));
		helper_simd.counters[TWK_LD_SIMD_ALTALT] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarA));
	}

	for(; k < this->byte_width; ++k){
		mask = ~(arrayA_mask[k] | arrayB_mask[k]);
		helper_simd.counters[TWK_LD_SIMD_REFREF] += POPCOUNT_ITER(((~arrayA[k]) & (~arrayB[k])) & mask);
		helper_simd.counters[TWK_LD_SIMD_REFALT] += POPCOUNT_ITER(((arrayA[k] ^ arrayB[k]) & arrayA[k]) & mask);
		helper_simd.counters[TWK_LD_SIMD_ALTREF] += POPCOUNT_ITER(((arrayA[k] ^ arrayB[k]) & arrayB[k]) & mask);
		helper_simd.counters[TWK_LD_SIMD_ALTALT] += POPCOUNT_ITER((arrayA[k] & arrayB[k]) & mask);
	}

	helper.alleleCounts[1] = helper_simd.counters[TWK_LD_SIMD_ALTREF];
	helper.alleleCounts[4] = helper_simd.counters[TWK_LD_SIMD_REFALT];
	helper.alleleCounts[5] = helper_simd.counters[TWK_LD_SIMD_ALTALT];
	helper.alleleCounts[0] = helper_simd.counters[TWK_LD_SIMD_REFREF] + (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 - phased_unbalanced_adjustment;
	//helper.haplotypeCounts[0] += (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT;


#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
	//std::cout << "V\t" << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif


	//this->setFLAGs(block1, block2);
	//return(this->CalculateLDPhasedMath());
	//std::cerr << "m1front=" << frontBonus << "," << tailBonus << "," << phased_unbalanced_adjustment << std::endl;
#if SLAVE_DEBUG_MODE == 2
	std::cerr << "m1 " << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[4] << "," << helper.alleleCounts[5] << std::endl;
#endif

#if SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::PhasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, twk_ld_perf* perf){
	helper.resetPhased();
	helper_simd.counters[TWK_LD_SIMD_ALTALT] = 0;

	const twk_igt_vec& block1 = b1.vec[p1];
	const twk_igt_vec& block2 = b2.vec[p2];
	const uint8_t* const arrayA = (const uint8_t* const)block1.data;
	const uint8_t* const arrayB = (const uint8_t* const)block2.data;

	// Debug timings
#if SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#if SIMD_AVAILABLE == 1
	const uint32_t frontSmallest = block1.front_zero < block2.front_zero ? block1.front_zero : block2.front_zero;
	const uint32_t tailSmallest  = block1.tail_zero  < block2.tail_zero  ? block1.tail_zero  : block2.tail_zero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	VECTOR_TYPE __intermediate;

#define ITER_SHORT {                                              \
	__intermediate  = PHASED_ALTALT(vectorA[i], vectorB[i]);      \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_ALTALT], __intermediate);\
	i += 1;                                                       \
}

	uint32_t i = frontSmallest;
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT


#undef ITER_SHORT
	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif

	//std::cerr << "k=" << k << "/" << this->byte_width << std::endl;
	for(; k+8 < this->byte_width; k += 8){
		for(uint32_t l = 0; l < 8; ++l)
			helper_simd.scalarA[l] = (arrayA[k+l] & arrayB[k+l]);
		helper_simd.counters[TWK_LD_SIMD_ALTALT] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarA));
	}

	for(; k < this->byte_width; ++k){
		helper_simd.counters[TWK_LD_SIMD_ALTALT] += POPCOUNT_ITER(arrayA[k] & arrayB[k]);
	}
	//helper.alleleCounts[5] = helper_simd.counters[TWK_LD_SIMD_REFREF] + (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 - this->phased_unbalanced_adjustment;

	helper.alleleCounts[5] = helper_simd.counters[TWK_LD_SIMD_ALTALT];
	helper.alleleCounts[4] = b1.blk->rcds[p1].ac - helper.alleleCounts[5];
	helper.alleleCounts[1] = b2.blk->rcds[p2].ac - helper.alleleCounts[5];
	helper.alleleCounts[0] = 2*n_samples - ((b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac) - helper.alleleCounts[5]);

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
	//std::cout << "V\t" << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << "\t" << ticks_per_iter.count() << '\n';
#endif

#if SLAVE_DEBUG_MODE == 2
	std::cerr << "m3=" << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[4] << "," << helper.alleleCounts[5] << std::endl;
#endif

	//this->setFLAGs(block1, block2);
	//return(this->CalculateLDPhasedMathSimple(block1, block2));
#if SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::UnphasedVectorized(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, twk_ld_perf* perf){
#if SLAVE_DEBUG_MODE == 0
	if(b1.blk->rcds[p1].gt_missing == false && b2.blk->rcds[p2].gt_missing == false){
		return(this->UnphasedVectorizedNoMissing(b1,p1,b2,p2,perf));
	}
#endif

	helper.resetUnphased();

	helper_simd.counters[0] = 0;
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;
	helper_simd.counters[4] = 0;
	helper_simd.counters[5] = 0;
	helper_simd.counters[6] = 0;
	helper_simd.counters[7] = 0;
	helper_simd.counters[8] = 0;
	helper_simd.counters[9] = 0;

	// Data
	const twk_igt_vec& block1 = b1.vec[p1];
	const twk_igt_vec& block2 = b2.vec[p2];
	const uint8_t* const arrayA = (const uint8_t* const)block1.data;
	const uint8_t* const arrayB = (const uint8_t* const)block2.data;
	const uint8_t* const arrayA_mask = (const uint8_t* const)block1.mask;
	const uint8_t* const arrayB_mask = (const uint8_t* const)block2.mask;

#if SIMD_AVAILABLE == 1
	const uint32_t frontSmallest = block1.front_zero < block2.front_zero ? block1.front_zero : block2.front_zero;
	const uint32_t tailSmallest  = block1.tail_zero  < block2.tail_zero  ? block1.tail_zero  : block2.tail_zero;
	const uint32_t frontBonus    = block1.front_zero != frontSmallest ? block1.front_zero : block2.front_zero;
	const uint32_t tailBonus     = block1.tail_zero  != tailSmallest  ? block1.tail_zero  : block2.tail_zero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	const VECTOR_TYPE* const vectorA_mask = (const VECTOR_TYPE* const)arrayA_mask;
	const VECTOR_TYPE* const vectorB_mask = (const VECTOR_TYPE* const)arrayB_mask;

	VECTOR_TYPE __intermediate, mask;

	uint32_t i = frontSmallest;
	VECTOR_TYPE altalt, refref, altref, refalt;
	//VECTOR_TYPE t00, t05, t50, t55, combTop, combLeft, combRight, combBottom;


// Debug timings
#if SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_BASE {                                             \
	mask	= MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);     \
	refref  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], mask);	\
	refalt  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], mask);	\
	altref  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], mask);	\
	altalt  = PHASED_ALTALT_MASK(vectorA[i], vectorB[i], mask);	\
	i += 1;                                                     \
}

#define ITER_SHORT {                                                        \
	ITER_BASE                                                               \
	__intermediate = FILTER_UNPHASED_SPECIAL(refref);                       \
	popcnt128(helper_simd.counters[0], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);    \
	popcnt128(helper_simd.counters[1], __intermediate);                      \
	__intermediate = FILTER_UNPHASED(altref, altref);                       \
	popcnt128(helper_simd.counters[2], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);    \
	popcnt128(helper_simd.counters[3], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refref, altalt, altalt, refref);  \
	popcnt128(helper_simd.counters[4], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altref, altref, refalt);  \
	popcnt128(helper_simd.counters[4], __intermediate);                      \
	__intermediate = FILTER_UNPHASED(refalt,refalt);                        \
	popcnt128(helper_simd.counters[6], __intermediate);                      \
}

#define ITER_LONG {                                                         \
	ITER_SHORT                                                              \
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref);    \
	popcnt128(helper_simd.counters[5], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt);  \
	popcnt128(helper_simd.counters[7], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt);                       \
	popcnt128(helper_simd.counters[8], __intermediate);                      \
}

	for( ; i < frontBonus; ) 					  	 ITER_SHORT
	for( ; i < this->vector_cycles - tailBonus; )  	 ITER_LONG
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT

#undef ITER_LONG
#undef ITER_SHORT
#undef ITER_BASE

	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif

	uint8_t b_altalt, b_refref, b_refalt, b_altref;
	for(; k < this->byte_width; ++k){
		b_altalt  = (arrayA[k] & arrayB[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_refref  = ((~arrayA[k]) & (~arrayB[k])) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_altref  = ((arrayA[k] ^ arrayB[k]) & arrayA[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_refalt  = ((arrayA[k] ^ arrayB[k]) & arrayB[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);

		helper_simd.counters[0] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_SPECIAL(b_refref));
		helper_simd.counters[1] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_refalt, b_refalt, b_refref));
		helper_simd.counters[2] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE(b_refalt, b_refalt));
		helper_simd.counters[3] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_altref, b_altref, b_refref));
		helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_altalt, b_altalt, b_refref));
		helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refalt, b_altref, b_altref, b_refalt));
		helper_simd.counters[5] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refalt, b_altalt, b_altalt, b_refalt));
		helper_simd.counters[6] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE(b_altref, b_altref));
		helper_simd.counters[7] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_altref, b_altalt, b_altalt, b_altref));
		helper_simd.counters[8] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_SPECIAL(b_altalt));
	}

	helper.alleleCounts[0]  = helper_simd.counters[0] - this->unphased_unbalanced_adjustment;
	helper.alleleCounts[0] += (frontSmallest + tailSmallest) * GENOTYPE_TRIP_COUNT;
	helper.alleleCounts[1]  = helper_simd.counters[1];
	helper.alleleCounts[5]  = helper_simd.counters[2];
	helper.alleleCounts[16] = helper_simd.counters[3];
	helper.alleleCounts[17] = helper_simd.counters[4];
	helper.alleleCounts[21] = helper_simd.counters[5];
	helper.alleleCounts[80] = helper_simd.counters[6];
	helper.alleleCounts[81] = helper_simd.counters[7];
	helper.alleleCounts[85] = helper_simd.counters[8];

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	//std::cout << ticks_per_iter.count() << '\n';
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
#endif

#if SLAVE_DEBUG_MODE == 2
	std::cerr << "miss1=" << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[5]
	          << "," << helper.alleleCounts[16] << "," << helper.alleleCounts[17] << "," << helper.alleleCounts[21]
	          << "," << helper.alleleCounts[80] << "," << helper.alleleCounts[81] << "," << helper.alleleCounts[85] << std::endl;
#endif
	//this->setFLAGs(block1, block2);
	//return(this->CalculateLDUnphasedMath());

#if SLAVE_DEBUG_MODE != 1
	return(UnphasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::UnphasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, twk_ld_perf* perf){
	helper.resetUnphased();

	helper_simd.counters[0] = 0;
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;
	helper_simd.counters[4] = 0;
	helper_simd.counters[5] = 0;
	helper_simd.counters[6] = 0;
	helper_simd.counters[7] = 0;
	helper_simd.counters[8] = 0;
	helper_simd.counters[9] = 0;

	// Data
	const twk_igt_vec& block1 = b1.vec[p1];
	const twk_igt_vec& block2 = b2.vec[p2];
	const uint8_t* const arrayA = (const uint8_t* const)block1.data;
	const uint8_t* const arrayB = (const uint8_t* const)block2.data;

#if SIMD_AVAILABLE == 1
	const uint32_t frontSmallest = block1.front_zero < block2.front_zero ? block1.front_zero : block2.front_zero;
	int i = frontSmallest;
	const uint32_t tailSmallest  = block1.tail_zero  < block2.tail_zero  ? block1.tail_zero  : block2.tail_zero;
	const uint32_t frontBonus    = block1.front_zero != frontSmallest ? block1.front_zero : block2.front_zero;
	const uint32_t tailBonus     = block1.tail_zero  != tailSmallest  ? block1.tail_zero  : block2.tail_zero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;

	VECTOR_TYPE altalt, refref, altref, refalt;
	VECTOR_TYPE __intermediate;

// Debug timings
#if SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_BASE {										\
	refref  = PHASED_REFREF(vectorA[i], vectorB[i]);	\
	refalt  = PHASED_REFALT(vectorA[i], vectorB[i]);	\
	altref  = PHASED_ALTREF(vectorA[i], vectorB[i]);	\
	altalt  = PHASED_ALTALT(vectorA[i], vectorB[i]);	\
	i += 1;												\
}

#define ITER_SHORT {														\
	ITER_BASE																\
	__intermediate = FILTER_UNPHASED_SPECIAL(refref);						\
	popcnt128(helper_simd.counters[0], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);	\
	popcnt128(helper_simd.counters[1], __intermediate);				\
	__intermediate = FILTER_UNPHASED(altref, altref);						\
	popcnt128(helper_simd.counters[2], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);	\
	popcnt128(helper_simd.counters[3], __intermediate);				\
	__intermediate = FILTER_UNPHASED(refalt,refalt);						\
	popcnt128(helper_simd.counters[6], __intermediate);				\
}

#define ITER_LONG {															\
	ITER_SHORT																\
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref);	\
	popcnt128(helper_simd.counters[5], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt);	\
	popcnt128(helper_simd.counters[7], __intermediate);				\
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt);						\
	popcnt128(helper_simd.counters[8], __intermediate);				\
}

	for( ; i < frontBonus; )                         ITER_SHORT
	for( ; i < this->vector_cycles - tailBonus; )    ITER_LONG
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT

#undef ITER_LONG
#undef ITER_SHORT
#undef ITER_BASE

	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif

	uint8_t b_altalt, b_refref, b_refalt, b_altref;
	for(; k < this->byte_width; ++k){
		b_altalt  = (arrayA[k] & arrayB[k]);
		b_refref  = ((~arrayA[k]) & (~arrayB[k]));
		b_altref  = ((arrayA[k] ^ arrayB[k]) & arrayA[k]);
		b_refalt  = ((arrayA[k] ^ arrayB[k]) & arrayB[k]);

		helper_simd.counters[0] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_SPECIAL(b_refref));
		helper_simd.counters[1] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_refalt, b_refalt, b_refref));
		helper_simd.counters[2] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE(b_refalt, b_refalt));
		helper_simd.counters[3] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_altref, b_altref, b_refref));
		//helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_altalt, b_altalt, b_refref));
		//helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refalt, b_altref, b_altref, b_refalt));

		helper_simd.counters[5] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refalt, b_altalt, b_altalt, b_refalt));
		helper_simd.counters[6] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE(b_altref, b_altref));
		helper_simd.counters[7] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_altref, b_altalt, b_altalt, b_altref));
		helper_simd.counters[8] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_SPECIAL(b_altalt));
	}

	helper.alleleCounts[0]  = helper_simd.counters[0] - this->unphased_unbalanced_adjustment;
	helper.alleleCounts[0] += (frontSmallest + tailSmallest) * GENOTYPE_TRIP_COUNT;
	helper.alleleCounts[1]  = helper_simd.counters[1];
	helper.alleleCounts[5]  = helper_simd.counters[2];
	helper.alleleCounts[16] = helper_simd.counters[3];
	helper.alleleCounts[21] = helper_simd.counters[5];
	helper.alleleCounts[80] = helper_simd.counters[6];
	helper.alleleCounts[81] = helper_simd.counters[7];
	helper.alleleCounts[85] = helper_simd.counters[8];
	//helper[17] = helper_simd.counters[4];
	helper.alleleCounts[17] = n_samples - (helper.alleleCounts[0] +  helper.alleleCounts[1] + helper.alleleCounts[5] + helper.alleleCounts[16] + helper.alleleCounts[21] + helper.alleleCounts[80] + helper.alleleCounts[81] + helper.alleleCounts[85]);

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	//std::cout << ticks_per_iter.count() << '\n';
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
#endif

#if SLAVE_DEBUG_MODE == 2
	std::cerr << "miss2=" << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[5]
		          << "," << helper.alleleCounts[16] << "," << helper.alleleCounts[17] << "," << helper.alleleCounts[21]
		          << "," << helper.alleleCounts[80] << "," << helper.alleleCounts[81] << "," << helper.alleleCounts[85] << std::endl;
#endif

	//this->setFLAGs(block1, block2);
#if SLAVE_DEBUG_MODE != 1
	return(UnphasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::PhasedRunlength(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, twk_ld_perf* perf){
	helper.resetPhased();
#if SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	const twk1_gt_t& gt1 = *b1.blk->rcds[p1].gt;
	const twk1_gt_t& gt2 = *b2.blk->rcds[p2].gt;

	uint32_t lenA = gt1.GetLength(0);
	uint32_t lenB = gt2.GetLength(0);
	uint8_t currentMixL = (gt1.GetRefA(0) << 2) | gt2.GetRefA(0);
	uint8_t currentMixR = (gt1.GetRefB(0) << 2) | gt2.GetRefB(0);
	uint32_t offsetA = 0;
	uint32_t offsetB = 0;
	uint32_t add;

	while(true){
		if(lenA > lenB){ // If processed run length A > processed run length B
			lenA -= lenB;
			add = lenB;
			lenB = gt2.GetLength(++offsetB);
		} else if(lenA < lenB){ // If processed run length A < processed run length B
			lenB -= lenA;
			add = lenA;
			lenA = gt1.GetLength(++offsetA);
		} else { // If processed run length A == processed run length B
			add = lenB;
			lenA = gt1.GetLength(++offsetA);
			lenB = gt2.GetLength(++offsetB);
		}
		helper.alleleCounts[currentMixL] += add;
		helper.alleleCounts[currentMixR] += add;

		// Exit condition
		if(offsetA == gt1.n || offsetB == gt2.n){
			if(offsetA != gt1.n || offsetB != gt2.n){
				std::cerr << utility::timestamp("FATAL") << "Failed to exit equally!\n" << offsetA << "/" << gt1.n << " and " << offsetB << "/" << gt2.n << std::endl;
				exit(1);
			}
			break;
		}

		currentMixL = (gt1.GetRefA(offsetA) << 2) | gt2.GetRefA(offsetB);
		currentMixR = (gt1.GetRefB(offsetA) << 2) | gt2.GetRefB(offsetB);
	}

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[gt1.n + gt2.n] += ticks_per_iter.count();
	++perf->freq[gt1.n + gt2.n];
	//std::cout << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif

#if SLAVE_DEBUG_MODE == 2
	std::cerr << "rleP=" << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[4] << "," << helper.alleleCounts[5] << std::endl;
#endif

	//std::cerr << "p= " << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[4] << "," << helper.alleleCounts[5] << std::endl;

#if SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::UnphasedRunlength(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, twk_ld_perf* perf){
	helper.resetUnphased();
#if SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	const twk1_gt_t& gt1 = *b1.blk->rcds[p1].gt;
	const twk1_gt_t& gt2 = *b2.blk->rcds[p2].gt;

	uint32_t lenA = gt1.GetLength(0);
	uint32_t lenB = gt2.GetLength(0);
	uint8_t  currentMix = (gt1.GetRefA(0) << 6) | (gt1.GetRefB(0) << 4) | (gt2.GetRefA(0) << 2) | (gt2.GetRefB(0));
	uint32_t offsetA = 0;
	uint32_t offsetB = 0;
	uint32_t add;

	while(true){
		if(lenA > lenB){ // If processed run length A > processed run length B
			lenA -= lenB;
			add = lenB;
			lenB = gt2.GetLength(++offsetB);
		} else if(lenA < lenB){ // If processed run length A < processed run length B
			lenB -= lenA;
			add = lenA;
			lenA = gt1.GetLength(++offsetA);
		} else { // If processed run length A == processed run length B
			add = lenB;
			lenA = gt1.GetLength(++offsetA);
			lenB = gt2.GetLength(++offsetB);
		}
		helper.alleleCounts[currentMix] += add;

		// Exit condition
		if(offsetA == gt1.n || offsetB == gt2.n){
			if(offsetA != gt1.n || offsetB != gt2.n){
				std::cerr << utility::timestamp("FATAL") << "Failed to exit equally!\n" << offsetA << "/" << gt1.n << " and " << offsetB << "/" << gt2.n << std::endl;
				exit(1);
			}
			break;
		}

		currentMix = (gt1.GetRefA(offsetA) << 6) | (gt1.GetRefB(offsetA) << 4) | (gt2.GetRefA(offsetB) << 2) | (gt2.GetRefB(offsetB));
	}

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[gt1.n + gt2.n] += ticks_per_iter.count();
	++perf->freq[gt1.n + gt2.n];
	//std::cout << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif

#if SLAVE_DEBUG_MODE == 2
	std::cerr << "rleU=" << helper.alleleCounts[0] << "," << helper.alleleCounts[1]+ helper.alleleCounts[4] << "," << helper.alleleCounts[5]
						 << "," << helper.alleleCounts[16]+helper.alleleCounts[64] << "," << helper.alleleCounts[17]+ helper.alleleCounts[20]+ helper.alleleCounts[65] + helper.alleleCounts[68] << "," << helper.alleleCounts[21]+ helper.alleleCounts[69]
						 << "," << helper.alleleCounts[80] << "," << helper.alleleCounts[81]+ helper.alleleCounts[84] << "," << helper.alleleCounts[85] << std::endl;
#endif

#if SLAVE_DEBUG_MODE != 1
	return(UnphasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::PhasedMath(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2){
	// Total amount of non-missing alleles
	helper.totalHaplotypeCounts = helper.alleleCounts[0] + helper.alleleCounts[1] +
								  helper.alleleCounts[4] + helper.alleleCounts[5];

	// All values are missing
	if(helper.totalHaplotypeCounts < MINIMUM_ALLOWED_ALLELES){
		//++this->insufficent_alleles;
		return false;
	}

	// Haplotype frequencies
	double pA = TWK_HAP_FREQ(helper,0); // pA
	double qA = TWK_HAP_FREQ(helper,1); // qA
	double pB = TWK_HAP_FREQ(helper,4); // pB
	double qB = TWK_HAP_FREQ(helper,5); // qB

	if(pA*qB - qA*pB == 0){ return false; }

	// Allelic frequencies
	const double g0 = ((double)helper.alleleCounts[0] + helper.alleleCounts[1]) / (helper.totalHaplotypeCounts);
	const double g1 = ((double)helper.alleleCounts[4] + helper.alleleCounts[5]) / (helper.totalHaplotypeCounts);
	const double h0 = ((double)helper.alleleCounts[0] + helper.alleleCounts[4]) / (helper.totalHaplotypeCounts);
	const double h1 = ((double)helper.alleleCounts[1] + helper.alleleCounts[5]) / (helper.totalHaplotypeCounts);

	cur_rcd2.D = pA*qB - qA*pB;
	cur_rcd2.R2 = cur_rcd2.D*cur_rcd2.D / (g0*g1*h0*h1);
	if(cur_rcd2.R2 < settings.minR2 || cur_rcd2.R2 > settings.maxR2){
		cur_rcd2.controller = 0;
		return false;
	}

	double dmax = 0;
	if(cur_rcd2.D >= 0) dmax = g0*h1 < h0*g1 ? g0*h1 : h0*g1;
	else dmax = g0*g1 < h0*h1 ? -g0*g1 : -h0*h1;

	cur_rcd2.Dprime = cur_rcd2.D / dmax;

	if(cur_rcd2.Dprime < settings.minDprime || cur_rcd2.Dprime > settings.maxDprime){
		cur_rcd2.controller = 0;
		return false;
	}

	// Calculate Fisher's exact test P-value
	double left,right,both;
	kt_fisher_exact(helper.alleleCounts[0],
					helper.alleleCounts[1],
					helper.alleleCounts[4],
					helper.alleleCounts[5],
					&left, &right, &both);

	if(both > settings.minP){
		cur_rcd2.controller = 0;
		return false;
	}
	cur_rcd2.P = both;
	cur_rcd2.R = sqrt(cur_rcd2.R2);

	cur_rcd2.Apos = b1.blk->rcds[p1].pos;
	cur_rcd2.Bpos = b2.blk->rcds[p2].pos;
	cur_rcd2.ridA = b1.blk->rcds[p1].rid;
	cur_rcd2.ridB = b2.blk->rcds[p2].rid;
	cur_rcd2[0] = helper.alleleCounts[0];
	cur_rcd2[1] = helper.alleleCounts[1];
	cur_rcd2[2] = helper.alleleCounts[4];
	cur_rcd2[3] = helper.alleleCounts[5];

	cur_rcd2.SetLowACA(b1.blk->rcds[p1].ac < LOW_AC_THRESHOLD);
	cur_rcd2.SetLowACB(b2.blk->rcds[p2].ac < LOW_AC_THRESHOLD);
	cur_rcd2.SetCompleteLD(helper.alleleCounts[0] < 1 || helper.alleleCounts[1] < 1 || helper.alleleCounts[4] < 1 || helper.alleleCounts[5] < 1);
	cur_rcd2.SetPerfectLD(cur_rcd2.R2 > 0.99);
	cur_rcd2.SetHasMissingValuesA(b1.blk->rcds[p1].an);
	cur_rcd2.SetHasMissingValuesB(b2.blk->rcds[p2].an);
	int32_t diff = (int32_t)b1.blk->rcds[p1].pos - b2.blk->rcds[p2].pos;
	cur_rcd2.SetLongRange(abs(diff) > LONG_RANGE_THRESHOLD && (cur_rcd2.ridA == cur_rcd2.ridB));
	cur_rcd2.SetUsedPhasedMath();
	cur_rcd2.SetSameContig(cur_rcd.ridA == cur_rcd.ridB);

	// Calculate Chi-Sq CV from 2x2 contingency table
	cur_rcd2.ChiSqModel = 0;
	cur_rcd2.ChiSqFisher = helper.totalHaplotypeCounts * cur_rcd2.R2;

#if SLAVE_DEBUG_MODE == 3
	if(cur_rcd2.R2 > 0.1){
		twk_debug_pos1 = b1.blk->rcds[p1].pos;
		twk_debug_pos1_2 = b2.blk->rcds[p2].pos;
//std::cout << "P\t" << b1.blk->rcds[p1].pos << "\t" << b2.blk->rcds[p2].pos << "\t" << helper.D << "\t" << helper.Dprime << "\t" << helper.R << "\t" << helper.R2 << '\n';
	}
#endif
	// If the number of rcds written is equal to the flush limit then
	// compress and write output.
	if(n_out == n_lim || irecF.rid != cur_rcd2.ridA || irecR.rid != cur_rcd2.ridB){
		if(this->CompressBlock() == false) return false;
		irecF.rid    = cur_rcd2.ridA;
		irecF.ridB   = cur_rcd2.ridB;
		irecF.minpos = cur_rcd2.Apos;
		irecF.maxpos = cur_rcd2.Apos;
		irecR.rid    = cur_rcd2.ridB;
		irecR.ridB   = cur_rcd2.ridA;
		irecR.minpos = cur_rcd2.Bpos;
		irecR.maxpos = cur_rcd2.Bpos;
	}
	// If ridB is mixed then set index value to -1.
	if(irecF.ridB != cur_rcd2.ridB) irecF.ridB = -1;
	if(irecR.ridB != cur_rcd2.ridA) irecR.ridB = -1;

	// Update index
	irecF.maxpos = cur_rcd2.Apos;
	irecR.maxpos = cur_rcd2.Bpos;

	// Add forward.
	blk_f += cur_rcd2;
	// Swap tuple (ridA:posA) with (ridB:posB).
	uint32_t temp = cur_rcd2.Apos;
	cur_rcd2.Apos = cur_rcd2.Bpos;
	cur_rcd2.Bpos = temp;
	std::swap(cur_rcd2.ridA,cur_rcd2.ridB);
	// Add reverse.
	blk_r += cur_rcd2;
	// Update tickers.
	n_out += 2; n_out_tick += 2;

	// Update ticker.
	if(n_out_tick == 300){ progress->n_out += 300; n_out_tick = 0; }

	// Reset controller.
	cur_rcd2.controller = 0;

	return true;
}

bool twk_ld_engine::UnphasedMath(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2){
	// Total amount of non-missing alleles
	helper.totalHaplotypeCounts = helper.alleleCounts[0]  + helper.alleleCounts[1]  + helper.alleleCounts[4]  + helper.alleleCounts[5]
	                     + helper.alleleCounts[16] + helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[21]
	                     + helper.alleleCounts[64] + helper.alleleCounts[65] + helper.alleleCounts[68] + helper.alleleCounts[69]
	                     + helper.alleleCounts[80] + helper.alleleCounts[81] + helper.alleleCounts[84] + helper.alleleCounts[85];

	// All values are missing or too few
	//if(helper.totalHaplotypeCounts < MINIMUM_ALLOWED_ALLELES){
		//++this->insufficent_alleles;
		//return false;
	//}

	// How many hets-hets is there? 0/1 0/1 or 1/0 1/0 or equivalent
	const uint64_t number_of_hets = helper.alleleCounts[17] + helper.alleleCounts[20]
	                              + helper.alleleCounts[65] + helper.alleleCounts[68];

	// If het-hets are 0 then do normal calculations
	// There is no phase uncertainty
	// Use phased math
#if SLAVE_DEBUG_MODE != 2
	if(number_of_hets == 0){
		const uint64_t p0 = 2*helper.alleleCounts[0] + helper.alleleCounts[1]  + helper.alleleCounts[4]    + helper.alleleCounts[16] + helper.alleleCounts[64];
		const uint64_t q0 = helper.alleleCounts[16]  + helper.alleleCounts[64] + 2*helper.alleleCounts[80] + helper.alleleCounts[81] + helper.alleleCounts[84];
		const uint64_t p1 = helper.alleleCounts[1]   + helper.alleleCounts[4]  + 2*helper.alleleCounts[5]  + helper.alleleCounts[21] + helper.alleleCounts[69];
		const uint64_t q1 = helper.alleleCounts[21]  + helper.alleleCounts[69] + helper.alleleCounts[81]   + helper.alleleCounts[84] + 2*helper.alleleCounts[85];

		helper.alleleCounts[0] = p0;
		helper.alleleCounts[1] = p1;
		helper.alleleCounts[4] = q0;
		helper.alleleCounts[5] = q1;

		// Update counter
		//++this->no_uncertainty;

		// Reset
		cur_rcd.ChiSqModel = 0;

		// Use standard math
		return(this->PhasedMath(b1,p1,b2,p2));
	}
#endif

	const double p = ((helper.alleleCounts[0] + helper.alleleCounts[1]  + helper.alleleCounts[4]  + helper.alleleCounts[5])*2.0
				   + (helper.alleleCounts[16] + helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[21] + helper.alleleCounts[64] + helper.alleleCounts[65] + helper.alleleCounts[68] + helper.alleleCounts[69]))
				   / (2.0 * helper.totalHaplotypeCounts);
	const double q = ((helper.alleleCounts[0] + helper.alleleCounts[16] + helper.alleleCounts[64] + helper.alleleCounts[80])*2.0
				   + (helper.alleleCounts[1]  + helper.alleleCounts[4]  + helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68] + helper.alleleCounts[81] + helper.alleleCounts[84]))
				   / (2.0 * helper.totalHaplotypeCounts);
	const double n11 = 2.0* helper.alleleCounts[0] + helper.alleleCounts[1] + helper.alleleCounts[4] + helper.alleleCounts[16] + helper.alleleCounts[64];

	// Not used for anything
	//const double n12 = (2.0*helper[5]  + helper[1]  + helper[4]  + helper[21] + helper[69]);
	//const double n21 = (2.0*helper[80] + helper[81] + helper[84] + helper[16] + helper[64]);
	//const double n22 = (2.0*helper[85] + helper[81] + helper[84] + helper[21] + helper[85]);

	/*////////////////////////
	// Cubic function: a3x^3 + a2x^2 + a1x + d = 0 <==> ax^3 + bx^2 + cx + d = 0
	// Cubic constants
	////////////////////////*/
	const double G   = 1.0 - 2.0*p - 2.0*q;
	const double dee = -n11*p*q;
	const double c   = -n11*G - number_of_hets*(1.0 - p - q) + 2.0*helper.totalHaplotypeCounts*p*q;
	const double b   = 2.0*helper.totalHaplotypeCounts*G - 2.0*n11 - number_of_hets;
	const double a   = 4.0 * helper.totalHaplotypeCounts;

	// Bounds for biological relevance
	const double minhap = n11 / (2.0 * helper.totalHaplotypeCounts);
	const double maxhap = (n11 + number_of_hets) / (2.0 * helper.totalHaplotypeCounts);

	// Cubic parameters
	const double xN  = -b / (3.0*a);
	const double d2  = (pow(b,2) - 3.0*a*c) / (9*pow(a,2));
	const double yN  = a * pow(xN,3) + b * pow(xN,2) + c * xN + dee;
	const double yN2 = pow(yN,2);
	const double h2  = 4 * pow(a,2) * pow(d2,3);

	// Difference between yN2 and h2
	const double __diff = yN2 - h2;

	// Begin cases
	if(__diff < 0) { // Yn2 < h2
		const double theta    = acos(-yN / sqrt(h2)) / 3.0;
		const double constant = 2.0 * sqrt(d2);
		const double alpha    = xN + constant * cos(theta);
		const double beta     = xN + constant * cos(2.0*M_PI/3.0 + theta);
		const double gamma    = xN + constant * cos(4.0*M_PI/3.0 + theta);

		uint8_t biologically_possible = 0;
		cur_rcd.ChiSqModel = std::numeric_limits<float>::max();
		const double* chosen = &alpha;
		if(alpha >= minhap - ALLOWED_ROUNDING_ERROR && alpha <= maxhap + ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(alpha, p, q);
			chosen = &alpha;
		}

		if(beta >= minhap - ALLOWED_ROUNDING_ERROR && beta <= maxhap + ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			if(this->ChiSquaredUnphasedTable(beta, p, q) < cur_rcd.ChiSqModel){
				chosen = &beta;
				cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(beta, p, q);
			}
		}

		if(gamma >= minhap - ALLOWED_ROUNDING_ERROR && gamma <= maxhap + ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			if(this->ChiSquaredUnphasedTable(gamma, p, q) < cur_rcd.ChiSqModel){
				chosen = &gamma;
				cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(gamma, p, q); // implicit
			}
		}

		if(biologically_possible == 0){
			//++this->impossible;
			return false;
		}

		if(biologically_possible > 1)
			cur_rcd.SetMultipleRoots();

		//++this->possible;
		return(this->ChooseF11Calculate(b1,p1,b2,p2,*chosen, p, q));

	} else if(__diff > 0){ // Yn2 > h2
		const double constant = sqrt(yN2 - h2);
		double left, right;
		if(1.0/(2.0*a)*(-yN + constant) < 0)
			 left = -pow(-1.0/(2.0*a)*(-yN + constant), 1.0/3.0);
		else left =  pow( 1.0/(2.0*a)*(-yN + constant), 1.0/3.0);
		if(1.0/(2.0*a)*(-yN - constant) < 0)
			 right = -pow(-1.0/(2.0*a)*(-yN - constant), 1.0/3.0);
		else right =  pow( 1.0/(2.0*a)*(-yN - constant), 1.0/3.0);

		const double alpha = xN + right + left;
		if(!(alpha >= minhap - ALLOWED_ROUNDING_ERROR && alpha <= maxhap + ALLOWED_ROUNDING_ERROR)){
			//++this->impossible;
			return false;
		}

		cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(alpha, p, q);
		//++this->possible;

		return(this->ChooseF11Calculate(b1,p1,b2,p2,alpha, p, q));

	} else { // Yn2 == h2
		const double delta = pow((yN/2.0*a),(1.0/3.0));
		const double alpha = xN + delta; // alpha = beta in this case
		const double gamma = xN - 2.0*delta;

		if(std::isnan(alpha) || std::isnan(gamma)){
			//++this->impossible;
			return false;
		}

		uint8_t biologically_possible = 0;
		cur_rcd.ChiSqModel = std::numeric_limits<float>::max();
		const double* chosen = &alpha;
		if(alpha >= minhap - ALLOWED_ROUNDING_ERROR && alpha <= maxhap + ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(alpha, p, q);
			chosen = &alpha;
		}

		if(gamma >= minhap - ALLOWED_ROUNDING_ERROR && gamma <= maxhap + ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			if(this->ChiSquaredUnphasedTable(gamma, p, q) < cur_rcd.ChiSqModel){
				chosen = &gamma;
				cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(gamma, p, q); // implicit
			}
		}

		if(biologically_possible == 0){
			//++this->impossible;
			return false;
		}

		//++this->possible;
		return(this->ChooseF11Calculate(b1,p1,b2,p2,*chosen, p, q));
	}
	return(false);
}

double twk_ld_engine::ChiSquaredUnphasedTable(const double& target, const double& p, const double& q){
	const double f12   = p - target;
	const double f21   = q - target;
	const double f22   = 1 - (target + f12 + f21);

	const double e1111 = helper.totalHaplotypeCounts * pow(target,2);
	const double e1112 = 2 * helper.totalHaplotypeCounts * target * f12;
	const double e1122 = helper.totalHaplotypeCounts * pow(f12,2);
	const double e1211 = 2 * helper.totalHaplotypeCounts * target * f21;
	const double e1212 = 2 * helper.totalHaplotypeCounts * f12 * f21 + 2 * helper.totalHaplotypeCounts * target * f22;
	const double e1222 = 2 * helper.totalHaplotypeCounts * f12 * f22;
	const double e2211 = helper.totalHaplotypeCounts * pow(f21,2);
	const double e2212 = 2 * helper.totalHaplotypeCounts * f21 * f22;
	const double e2222 = helper.totalHaplotypeCounts * pow(f22,2);

	const double chisq1111 = e1111 > 0 ? pow((double)helper.alleleCounts[0]  - e1111, 2) / e1111 : 0,
				 chisq1112 = e1112 > 0 ? pow((double)helper.alleleCounts[1]  + helper.alleleCounts[4] - e1112, 2) / e1112 : 0,
				 chisq1122 = e1122 > 0 ? pow((double)helper.alleleCounts[5]  - e1122, 2) / e1122 : 0,
				 chisq1211 = e1211 > 0 ? pow((double)helper.alleleCounts[16] + helper.alleleCounts[64] - e1211, 2) / e1211 : 0,
				 chisq1212 = e1212 > 0 ? pow((double)helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68] - e1212, 2) / e1212 : 0,
				 chisq1222 = e1222 > 0 ? pow((double)helper.alleleCounts[21] + helper.alleleCounts[69] - e1222, 2) / e1222 : 0,
				 chisq2211 = e2211 > 0 ? pow((double)helper.alleleCounts[80] - e2211, 2) / e2211 : 0,
				 chisq2212 = e2212 > 0 ? pow((double)helper.alleleCounts[81] + helper.alleleCounts[84] - e2212, 2) / e2212 : 0,
				 chisq2222 = e2222 > 0 ? pow((double)helper.alleleCounts[85] - e2222, 2) / e2222 : 0;

	return(chisq1111 + chisq1112 + chisq1122 + chisq1211 + chisq1212 + chisq1222 + chisq2211 + chisq2212 + chisq2222);
}

bool twk_ld_engine::ChooseF11Calculate(const twk1_ldd_blk& b1, const uint32_t& pos1, const twk1_ldd_blk& b2, const uint32_t& pos2,const double& target, const double& p, const double& q){
	const double p1 = target;
	const double p2 = p - p1;
	const double q1 = q - p1;
	const double q2 = (1-(p1+p2+q1) < 0 ? 0 : 1-(p1+p2+q1));

	//std::cerr << "unphased-choose=" << p1 << "," << p2 << "," << q1 << "," << q2 << " with p=" << p << " q=" << q  << std::endl;

	cur_rcd.D  = p1*q2 - p2*q1;
	cur_rcd.R2 = cur_rcd.D*cur_rcd.D / (p*(1-p)*q*(1-q));

	if(cur_rcd.R2 < settings.minR2 || cur_rcd.R2 > settings.maxR2){
		cur_rcd.controller = 0;
		return false;
	}

	cur_rcd.R  = sqrt(cur_rcd.R2);

	cur_rcd[0] = p1 * 2*helper.totalHaplotypeCounts;
	cur_rcd[1] = p2 * 2*helper.totalHaplotypeCounts;
	cur_rcd[2] = q1 * 2*helper.totalHaplotypeCounts;
	cur_rcd[3] = q2 * 2*helper.totalHaplotypeCounts;

	//std::cerr << p1 << "," << p2 << "," << q1 << "," << q2 << " and " << cur_rcd[0] << "," << cur_rcd[1] << "," << cur_rcd[2] << "," << cur_rcd[3] << " with R2=" << cur_rcd.R2 << std::endl;

	if(cur_rcd[1] + cur_rcd[2] + cur_rcd[3] <= 2){
		//std::cerr << "filter out=" << cur_rcd[1] + cur_rcd[2] + cur_rcd[3] << " when " << b1[pos1].ac << "," << b2[pos2].ac << std::endl;
		cur_rcd.controller = 0;
		return false;
	}

	double dmax = 0;
	if(cur_rcd.D >= 0) dmax = p*(1.0-q) < q*(1.0-p) ? p*(1.0-q) : q*(1.0-p);
	else dmax = p*q < (1-p)*(1-q) ? -p*q : -(1-p)*(1-q);
	cur_rcd.Dprime = cur_rcd.D / dmax;

	if(cur_rcd.Dprime < settings.minDprime || cur_rcd.Dprime > settings.maxDprime){
		cur_rcd.controller = 0;
		return false;
	}

	double left,right,both;
	kt_fisher_exact(round(cur_rcd[0]),round(cur_rcd[1]),
					round(cur_rcd[2]),round(cur_rcd[3]),
					&left,&right,&both);
	cur_rcd.P = both;

	if(cur_rcd.P > settings.minP){
		cur_rcd.controller = 0;
		return false;
	}

	cur_rcd.Apos = b1.blk->rcds[pos1].pos;
	cur_rcd.Bpos = b2.blk->rcds[pos2].pos;
	cur_rcd.ridA = b1.blk->rcds[pos1].rid;
	cur_rcd.ridB = b2.blk->rcds[pos2].rid;
	cur_rcd.ChiSqModel = 0;
	cur_rcd.ChiSqFisher = (cur_rcd[0] + cur_rcd[1] + cur_rcd[2] + cur_rcd[3]) * cur_rcd.R2;

	//std::cerr << "unphased -> D=" << helper.D << ",Dprime=" << helper.Dprime << ",R2=" << helper.R2 << ",R=" << helper.R << std::endl;
	cur_rcd.SetLowACA(b1.blk->rcds[pos1].ac < LOW_AC_THRESHOLD);
	cur_rcd.SetLowACB(b2.blk->rcds[pos2].ac < LOW_AC_THRESHOLD);
	cur_rcd.SetCompleteLD(cur_rcd[0] < 1 || cur_rcd[1] < 1 || cur_rcd[4] < 1 || cur_rcd[5] < 1);
	cur_rcd.SetPerfectLD(cur_rcd2.R2 > 0.99);
	cur_rcd.SetHasMissingValuesA(b1.blk->rcds[pos1].an);
	cur_rcd.SetHasMissingValuesB(b2.blk->rcds[pos2].an);
	int32_t diff = (int32_t)b1.blk->rcds[pos1].pos - b2.blk->rcds[pos2].pos;
	cur_rcd.SetLongRange(abs(diff) > LONG_RANGE_THRESHOLD && (cur_rcd.ridA == cur_rcd.ridB));
	cur_rcd.SetUsedPhasedMath();
	cur_rcd.SetSameContig(cur_rcd.ridA == cur_rcd.ridB);

	//if(helper.R2 > 0.5) exit(1);

#if SLAVE_DEBUG_MODE == 3
	if(cur_rcd.R2 > 0.1){
		twk_debug_pos2   = b1.blk->rcds[pos1].pos;
		twk_debug_pos2_2 = b2.blk->rcds[pos2].pos;
		if(twk_debug_pos1 == twk_debug_pos2 && twk_debug_pos1_2 == twk_debug_pos2_2){
			std::cout << "P\t" << cur_rcd2 << '\n';
			std::cout << "U\t" << cur_rcd << '\n';
		}
		//std::cout << "U\t" << b1.blk->rcds[pos1].pos << "\t" << b2.blk->rcds[pos2].pos << "\t" << helper.D << "\t" << helper.Dprime << "\t" << helper.R << "\t" << helper.R2 << '\n';
	}
	#endif

	// If the number of rcds written is equal to the flush limit then
	// compress and write output.
	if(n_out == n_lim || irecF.rid != cur_rcd.ridA || irecR.rid != cur_rcd.ridB){
		if(this->CompressBlock() == false) return false;
		irecF.rid    = cur_rcd.ridA;
		irecF.ridB   = cur_rcd.ridB;
		irecF.minpos = cur_rcd.Apos;
		irecF.maxpos = cur_rcd.Apos;
		irecR.rid    = cur_rcd.ridB;
		irecR.ridB   = cur_rcd.ridA;
		irecR.minpos = cur_rcd.Bpos;
		irecR.maxpos = cur_rcd.Bpos;
	}
	// If ridB is mixed then set index value to -1.
	if(irecF.ridB != cur_rcd.ridB) irecF.ridB = -1;
	if(irecR.ridB != cur_rcd.ridA) irecR.ridB = -1;

	// Update index
	irecF.maxpos = cur_rcd.Apos;
	irecR.maxpos = cur_rcd.Bpos;

	// Add forward.
	blk_f += cur_rcd;
	// Swap tuple (ridA:posA) with (ridB:posB).
	uint32_t temp = cur_rcd.Apos;
	cur_rcd.Apos = cur_rcd.Bpos;
	cur_rcd.Bpos = temp;
	std::swap(cur_rcd.ridA,cur_rcd.ridB);
	// Add reverse.
	blk_r += cur_rcd;
	// Update tickers.
	n_out += 2; n_out_tick += 2;

	// Update ticker.
	if(n_out_tick == 300){ progress->n_out += 300; n_out_tick = 0; }

	// Reset controller.
	cur_rcd.controller = 0;

	return false;
}

bool twk_ld_engine::CompressFwd(){
	if(blk_f.n){
		ibuf << blk_f;

		if(zcodec.Compress(ibuf, obuf, settings.c_level) == false){
			std::cerr << "failed compression" << std::endl;
			return false;
		}

		t_out += ibuf.size() / sizeof(twk1_two_t);
		progress->b_out += ibuf.size();
		//std::cerr << "F=" << buf_f.size() << "->" << obuf.size() << " -> " << (float)buf_f.size()/obuf.size() << std::endl;
		irecF.foff = writer->stream.tellp();
		writer->Add(ibuf.size(), obuf.size(), obuf);
		irecF.fend = writer->stream.tellp();
		irecF.n = blk_f.n;
		index->AddThreadSafe(irecF);
		//std::cerr << irecF.n << "," << irecF.foff << "," << irecF.fend << "," << irecF.rid << ":" << irecF.minpos << "-" << irecF.maxpos << "," << irecF.ridB << std::endl;
		ibuf.reset();
		blk_f.reset();
		irecF.clear();
	}

	return true;
}

bool twk_ld_engine::CompressRev(){
	if(blk_r.n){
		ibuf << blk_r;

		if(zcodec.Compress(ibuf, obuf, settings.c_level) == false){
			std::cerr << "failed compression" << std::endl;
			return false;
		}

		t_out += ibuf.size() / sizeof(twk1_two_t);
		progress->b_out += ibuf.size();
		//std::cerr << "F=" << buf_f.size() << "->" << obuf.size() << " -> " << (float)buf_f.size()/obuf.size() << std::endl;
		irecR.foff = writer->stream.tellp();
		writer->Add(ibuf.size(), obuf.size(), obuf);
		irecR.fend = writer->stream.tellp();
		irecR.n = blk_r.n;
		index->AddThreadSafe(irecR);
		//std::cerr << irecF.n << "," << irecF.foff << "," << irecF.fend << "," << irecF.rid << ":" << irecF.minpos << "-" << irecF.maxpos << "," << irecF.ridB << std::endl;
		ibuf.reset();
		blk_r.reset();
		irecR.clear();
	}

	return true;
}

bool twk_ld_engine::CompressBlock(){
	if(CompressFwd() == false) return false;
	if(CompressRev() == false) return false;
	n_out = 0;

	return true;
}

}
