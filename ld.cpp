#include <chrono>

#include "ld.h"

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

bool LDEngine::PhasedSimple(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, ld_perf* perf){
	helper.resetPhased();
	const twk_igt_list& ref = b1.list[p1];
	const twk_igt_list& tgt = b2.list[p2];
	const uint32_t n_cycles = ref.l_list < tgt.l_list ? ref.l_list : tgt.l_list;
	const uint32_t n_total  = ref.l_list + tgt.l_list;
	uint32_t n_same = 0;

	// Debug timings
	#if SLAVE_DEBUG_MODE == 1
		typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
		auto t0 = std::chrono::high_resolution_clock::now();
	#endif

	if(ref.l_list >= tgt.l_list){
		for(uint32_t i = 0; i < n_cycles; ++i) n_same += ref.get(tgt.list[i]);
	} else {
		for(uint32_t i = 0; i < n_cycles; ++i) n_same += tgt.get(ref.list[i]);
	}

	helper.haplotypeCounts[0] = 2*n_samples - (n_total - n_same);
	//std::cerr << 2*this->n_samples - (n_total - n_same) << '\n';
	//assert(2*n_samples - (n_total - n_same) < 2*n_samples);
	//std::cerr << helper.haplotypeCounts[0] << " ";

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[n_cycles] += ticks_per_iter.count();
	++perf->freq[n_cycles];
	//std::cout << "V\t" << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << "\t" << ticks_per_iter.count() << '\n';
#endif

	//return(this->CalculateLDPhasedMathSimple(block1, block2));
	return(MathPhasedSimple(b1,p1,b2,p2));
}

bool LDEngine::PhasedVectorized(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, ld_perf* perf){
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
	POPCOUNT(helper_simd.counters[TWK_LD_REFREF], __intermediate);       \
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks); \
	POPCOUNT(helper_simd.counters[TWK_LD_ALTREF], __intermediate);       \
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks); \
	POPCOUNT(helper_simd.counters[TWK_LD_REFALT], __intermediate);       \
	i += 1;                                                              \
}

#define ITER {                                                           \
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);              \
	__intermediate  = PHASED_ALTALT_MASK(vectorA[i], vectorB[i], masks); \
	POPCOUNT(helper_simd.counters[TWK_LD_ALTALT], __intermediate);       \
	__intermediate  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], masks); \
	POPCOUNT(helper_simd.counters[TWK_LD_REFREF], __intermediate);       \
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks); \
	POPCOUNT(helper_simd.counters[TWK_LD_ALTREF], __intermediate);       \
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks); \
	POPCOUNT(helper_simd.counters[TWK_LD_REFALT], __intermediate);       \
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
		helper_simd.counters[TWK_LD_REFREF] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarB));
		helper_simd.counters[TWK_LD_REFALT] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarC));
		helper_simd.counters[TWK_LD_ALTREF] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarD));
		helper_simd.counters[TWK_LD_ALTALT] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarA));
	}

	for(; k < this->byte_width; ++k){
		mask = ~(arrayA_mask[k] | arrayB_mask[k]);
		helper_simd.counters[TWK_LD_REFREF] += POPCOUNT_ITER(((~arrayA[k]) & (~arrayB[k])) & mask);
		helper_simd.counters[TWK_LD_REFALT] += POPCOUNT_ITER(((arrayA[k] ^ arrayB[k]) & arrayA[k]) & mask);
		helper_simd.counters[TWK_LD_ALTREF] += POPCOUNT_ITER(((arrayA[k] ^ arrayB[k]) & arrayB[k]) & mask);
		helper_simd.counters[TWK_LD_ALTALT] += POPCOUNT_ITER((arrayA[k] & arrayB[k]) & mask);
	}

	helper.haplotypeCounts[TWK_LD_ALTREF] = helper_simd.counters[TWK_LD_ALTREF];
	helper.haplotypeCounts[TWK_LD_REFALT] = helper_simd.counters[TWK_LD_REFALT];
	helper.haplotypeCounts[TWK_LD_ALTALT] = helper_simd.counters[TWK_LD_ALTALT];
	helper.haplotypeCounts[TWK_LD_REFREF] = helper_simd.counters[TWK_LD_REFREF] + (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 - phased_unbalanced_adjustment;
	//helper.haplotypeCounts[0] += (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT;
	//std::cerr << "m1=" << tailSmallest << "," << frontSmallest << "," << helper_simd.counters[0] << "," << phased_unbalanced_adjustment << std::endl;


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
	//std::cerr << "m1 " << helper.haplotypeCounts[0] << "," << helper.haplotypeCounts[1] << "," << helper.haplotypeCounts[2] << "," << helper.haplotypeCounts[3] << std::endl;

	return(MathPhasedSimple(b1,p1,b2,p2));
}

bool LDEngine::PhasedVectorizedNoMissingNoTable(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, ld_perf* perf){
	helper.resetPhased();
	helper_simd.counters[TWK_LD_REFREF] = 0;

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
	__intermediate  = PHASED_REFREF(vectorA[i], vectorB[i]);      \
	POPCOUNT(helper_simd.counters[TWK_LD_REFREF], __intermediate);\
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
			helper_simd.scalarB[l] = ((~arrayA[k+l]) & (~arrayB[k+l]));
		helper_simd.counters[TWK_LD_REFREF] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarB));
	}

	for(; k < this->byte_width; ++k){
		helper_simd.counters[TWK_LD_REFREF] += POPCOUNT_ITER(((~arrayA[k]) & (~arrayB[k])) & 255);
	}
	helper.haplotypeCounts[TWK_LD_REFREF] = helper_simd.counters[TWK_LD_REFREF] + (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 - this->phased_unbalanced_adjustment;
	//std::cerr << "m3=" << helper.haplotypeCounts[0] << "," << helper.haplotypeCounts[1] << "," << helper.haplotypeCounts[2] << "," << helper.haplotypeCounts[3] << std::endl;

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
	//std::cout << "V\t" << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << "\t" << ticks_per_iter.count() << '\n';
#endif

	//this->setFLAGs(block1, block2);
	//return(this->CalculateLDPhasedMathSimple(block1, block2));
	return(MathPhasedSimple(b1,p1,b2,p2));
}

bool LDEngine::PhasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, ld_perf* perf){
	helper.resetPhased();
	helper_simd.counters[0] = 0;
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;

	const twk_igt_vec& block1 = b1.vec[p1];
	const twk_igt_vec& block2 = b2.vec[p2];
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

	uint32_t i = frontSmallest;
	//VECTOR_TYPE altalt, altref, refalt;

	// Debug timings
#if SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_SHORT {											\
	__intermediate  = PHASED_REFALT(vectorA[i], vectorB[i]);	\
	POPCOUNT(helper_simd.counters[TWK_LD_REFALT], __intermediate);	\
	__intermediate  = PHASED_ALTREF(vectorA[i], vectorB[i]);	\
	POPCOUNT(helper_simd.counters[TWK_LD_ALTREF], __intermediate);	\
	i += 1;														\
}

#define ITER {													\
	__intermediate  = PHASED_ALTALT(vectorA[i], vectorB[i]);	\
	POPCOUNT(helper_simd.counters[TWK_LD_ALTALT], __intermediate);	\
	ITER_SHORT													\
}

	for( ; i < frontBonus; ) 					  	ITER_SHORT
	for( ; i < this->vector_cycles - tailBonus; )  	ITER
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT

#undef ITER
#undef ITER_SHORT
	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif


#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k+8 < this->byte_width; k += 8){
#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
		for(uint32_t l = 0; l < 8; ++l){
			helper_simd.scalarA[l] = (arrayA[k+l] & arrayB[k+l]);
			helper_simd.scalarC[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayA[k+l]);
			helper_simd.scalarD[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayB[k+l]);
		}

		helper_simd.counters[TWK_LD_REFALT] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarC));
		helper_simd.counters[TWK_LD_ALTREF] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarD));
		helper_simd.counters[TWK_LD_ALTALT] += POPCOUNT_ITER(*reinterpret_cast<const uint64_t* const>(helper_simd.scalarA));
	}

#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k < this->byte_width; ++k){
		helper_simd.counters[TWK_LD_REFALT] += POPCOUNT_ITER((arrayA[k] ^ arrayB[k]) & arrayA[k]);
		helper_simd.counters[TWK_LD_ALTREF] += POPCOUNT_ITER((arrayA[k] ^ arrayB[k]) & arrayB[k]);
		helper_simd.counters[TWK_LD_ALTALT] += POPCOUNT_ITER(arrayA[k] & arrayB[k]);
	}

	helper.haplotypeCounts[TWK_LD_ALTREF] = helper_simd.counters[TWK_LD_ALTREF];
	helper.haplotypeCounts[TWK_LD_REFALT] = helper_simd.counters[TWK_LD_REFALT];
	helper.haplotypeCounts[TWK_LD_ALTALT] = helper_simd.counters[TWK_LD_ALTALT];
	helper.haplotypeCounts[TWK_LD_REFREF] = 2*n_samples - (helper.haplotypeCounts[1] + helper.haplotypeCounts[2] + helper.haplotypeCounts[3]);
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
	//std::cout << "V\t" << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << "\t" << ticks_per_iter.count() << '\n';
#endif

	//this->setFLAGs(block1, block2);
	//std::cerr << "m2 " << helper.haplotypeCounts[0] << "," << helper.haplotypeCounts[1] << "," << helper.haplotypeCounts[2] << "," << helper.haplotypeCounts[3] << std::endl;

	//return(this->CalculateLDPhasedMath());
	return(MathPhasedSimple(b1,p1,b2,p2));
}

bool LDEngine::UnphasedVectorized(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, ld_perf* perf){
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
	POPCOUNT(helper_simd.counters[0], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);    \
	POPCOUNT(helper_simd.counters[1], __intermediate);                      \
	__intermediate = FILTER_UNPHASED(altref, altref);                       \
	POPCOUNT(helper_simd.counters[2], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);    \
	POPCOUNT(helper_simd.counters[3], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refref, altalt, altalt, refref);  \
	POPCOUNT(helper_simd.counters[4], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altref, altref, refalt);  \
	POPCOUNT(helper_simd.counters[4], __intermediate);                      \
	__intermediate = FILTER_UNPHASED(refalt,refalt);                        \
	POPCOUNT(helper_simd.counters[6], __intermediate);                      \
}

#define ITER_LONG {                                                         \
	ITER_SHORT                                                              \
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref);    \
	POPCOUNT(helper_simd.counters[5], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt);  \
	POPCOUNT(helper_simd.counters[7], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt);                       \
	POPCOUNT(helper_simd.counters[8], __intermediate);                      \
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

	//std::cerr << "miss1=" << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[5]
	//          << "," << helper.alleleCounts[16] << "," << helper.alleleCounts[17] << "," << helper.alleleCounts[21]
	//          << "," << helper.alleleCounts[80] << "," << helper.alleleCounts[81] << "," << helper.alleleCounts[85] << std::endl;

	//this->setFLAGs(block1, block2);
	//return(this->CalculateLDUnphasedMath());

	return true;
}

bool LDEngine::UnphasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, ld_perf* perf){
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
	POPCOUNT(helper_simd.counters[0], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);	\
	POPCOUNT(helper_simd.counters[1], __intermediate);				\
	__intermediate = FILTER_UNPHASED(altref, altref);						\
	POPCOUNT(helper_simd.counters[2], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);	\
	POPCOUNT(helper_simd.counters[3], __intermediate);				\
	__intermediate = FILTER_UNPHASED(refalt,refalt);						\
	POPCOUNT(helper_simd.counters[6], __intermediate);				\
}

#define ITER_LONG {															\
	ITER_SHORT																\
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref);	\
	POPCOUNT(helper_simd.counters[5], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt);	\
	POPCOUNT(helper_simd.counters[7], __intermediate);				\
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt);						\
	POPCOUNT(helper_simd.counters[8], __intermediate);				\
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

	//std::cerr << "miss2=" << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[5]
	//	          << "," << helper.alleleCounts[16] << "," << helper.alleleCounts[17] << "," << helper.alleleCounts[21]
	//	          << "," << helper.alleleCounts[80] << "," << helper.alleleCounts[81] << "," << helper.alleleCounts[85] << std::endl;


	//this->setFLAGs(block1, block2);
	//return(this->CalculateLDUnphasedMath());
	return true;
}

bool LDEngine::Runlength(const twk1_ldd_blk& b1, const uint32_t& p1, const twk1_ldd_blk& b2, const uint32_t& p2, ld_perf* perf){
	helper.resetPhased();
#if SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	const twk1_gt_t& gt1 = *b1.blk->rcds[p1].gt;
	const twk1_gt_t& gt2 = *b2.blk->rcds[p2].gt;

	uint32_t currentLengthA = gt1.GetLength(0);
	uint32_t currentLengthB = gt2.GetLength(0);
	uint8_t currentMixL = (gt1.GetRefA(0) << 2) | gt2.GetRefA(0);
	uint8_t currentMixR = (gt1.GetRefB(0) << 2) | gt2.GetRefB(0);
	uint32_t pointerA = 0;
	uint32_t pointerB = 0;
	uint32_t add;

#if SLAVE_DEBUG_MODE == 3
	U64 iterations = 0;
#endif

	while(true){
		if(currentLengthA > currentLengthB){ // If processed run length A > processed run length B
			currentLengthA -= currentLengthB;
			add = currentLengthB;
			currentLengthB = gt2.GetLength(++pointerB);
		} else if(currentLengthA < currentLengthB){ // If processed run length A < processed run length B
			currentLengthB -= currentLengthA;
			add = currentLengthA;
			currentLengthA = gt1.GetLength(++pointerA);
		} else { // If processed run length A == processed run length B
			add = currentLengthB;
			currentLengthA = gt1.GetLength(++pointerA);
			currentLengthB = gt2.GetLength(++pointerB);
		}
		helper.alleleCounts[currentMixL] += add;
		helper.alleleCounts[currentMixR] += add;

		// Exit condition
		if(pointerA == gt1.n || pointerB == gt2.n){
			if(pointerA != gt1.n || pointerB != gt2.n){
				std::cerr << utility::timestamp("FATAL") << "Failed to exit equally!\n" << pointerA << "/" << gt1.n << " and " << pointerB << "/" << gt2.n << std::endl;
				exit(1);
			}
			break;
		}

		currentMixL = (gt1.GetRefA(pointerA) << 2) | gt2.GetRefA(pointerB);
		currentMixR = (gt1.GetRefB(pointerA) << 2) | gt2.GetRefB(pointerB);

	#if SLAVE_DEBUG_MODE == 3
		++iterations;
		std::cout << a.getMeta().runs << '\t' << b.getMeta().runs << '\t' << iterations << std::endl;
	#endif
	}

#if SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[gt1.n + gt2.n] += ticks_per_iter.count();
	++perf->freq[gt1.n + gt2.n];
	//std::cout << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif

	//std::cerr << "p= " << helper.alleleCounts[0] << "," << helper.alleleCounts[1] << "," << helper.alleleCounts[4] << "," << helper.alleleCounts[5] << std::endl;


	return(MathPhasedSimple(b1,p1,b2,p2));
}

}
