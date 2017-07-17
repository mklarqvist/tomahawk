#ifndef TOMAHAWK_TOMAHAWKCALCULATESLAVE_H_
#define TOMAHAWK_TOMAHAWKCALCULATESLAVE_H_

#include <cmath>
#include <thread>
#include <cassert>
#include <cmath>

#include "../algorithm/spinlock.h"
#include "TomahawkOutputLD.h"
#include "TomahawkCalculationWriter.h"
#include "../interface/ProgressBar.h"
#include "TomahawkCalcParameters.h"
#include "TomahawkReaderControllerManager.h"
#include "TomahawkCalculateSlaveOutputManager.h"
#include "../algorithm/FisherTest.h"
#include "../algorithm/LoadBalancerBlock.h"
#include "../algorithm/GenotypeBitPacker.h"
#include "../algorithm/TomahawkSlaveSIMDHelper.h"

// Method 1: None: Input-specified (default)
// Method 2: Phased Vectorized No-Missing
// Method 3:
// Method 4: Unphased regular and unphased vectorized
// Method 5:
// Method 6: All algorithms comparison
#define __DEBUG_MODE__	1

namespace Tomahawk{

// Parameter flags
#define LOW_MAF_THRESHOLD		0.01
#define LOW_HWE_THRESHOLD		1e-6
#define LONG_RANGE_THRESHOLD	500e3
#define MINIMUM_ALLOWED_ALLELES	5
#define CHI_SQ_MAX_CV			1300

// SIMD trigger
#if SIMD_AVAILABLE == 1

#define POPCOUNT_ITER	_popcnt64

#define UNPHASED_UPPER_MASK	170  // 10101010b
#define UNPHASED_LOWER_MASK	85   // 01010101b
#define FILTER_UNPHASED_BYTE(A, B)	(((((A & UNPHASED_UPPER_MASK) | (B & UNPHASED_LOWER_MASK)) & UNPHASED_LOWER_MASK) << 1) & A)
#define FILTER_UNPHASED_BYTE_PAIR(A, B, C, D)	((FILTER_UNPHASED_BYTE(A, B) >> 1) | FILTER_UNPHASED_BYTE(C, D))
#define FILTER_UNPHASED_BYTE_SPECIAL(A)			(((A >> 1) & A) & UNPHASED_LOWER_MASK)

#if SIMD_VERSION == 6 // AVX-512: UNTESTED
const VECTOR_TYPE ONE_MASK = _mm512_set1_epi8(255); // 11111111b
const VECTOR_TYPE maskUnphasedHigh = _mm512_set1_epi8(UNPHASED_UPPER_MASK);	// 10101010b
const VECTOR_TYPE maskUnphasedLow  = _mm512_set1_epi8(UNPHASED_LOWER_MASK);	// 01010101b

#define PHASED_ALTALT(A,B)	_mm512_and_si512(A, B)
#define PHASED_REFREF(A,B)	_mm512_and_si512(_mm512_xor_si512(A, ONE_MASK), _mm512_xor_si512(B, ONE_MASK))
#define PHASED_ALTREF(A,B)	_mm512_and_si512(_mm512_xor_si512(A, B), B)
#define PHASED_REFALT(A,B)	_mm512_and_si512(_mm512_xor_si512(A, B), A)
#define PHASED_ALTALT_MASK(A,B,M)	_mm512_and_si512(PHASED_ALTALT(A, B), M)
#define PHASED_REFREF_MASK(A,B,M)	_mm512_and_si512(PHASED_REFREF(A, B), M)
#define PHASED_ALTREF_MASK(A,B,M)	_mm512_and_si512(PHASED_ALTREF(A, B), M)
#define PHASED_REFALT_MASK(A,B,M)	_mm512_and_si512(PHASED_REFALT(A, B), M)
#define MASK_MERGE(A,B)		_mm512_xor_si512(_mm512_or_si512(A, B), ONE_MASK)

// Software intrinsic popcount
#define POPCOUNT(A, B) {										\
	__m256i tempA = _mm512_extracti64x4_epi64(B, 0);			\
	A += __builtin_popcountll(_mm256_extract_epi64(tempA, 0));	\
	A += __builtin_popcountll(_mm256_extract_epi64(tempA, 1)); 	\
	A += __builtin_popcountll(_mm256_extract_epi64(tempA, 2)); 	\
	A += __builtin_popcountll(_mm256_extract_epi64(tempA, 3)); 	\
	tempA = _mm512_extracti64x4_epi64(B, 1); 					\
	A += __builtin_popcountll(_mm256_extract_epi64(tempA, 0));	\
	A += __builtin_popcountll(_mm256_extract_epi64(tempA, 1)); 	\
	A += __builtin_popcountll(_mm256_extract_epi64(tempA, 2)); 	\
	A += __builtin_popcountll(_mm256_extract_epi64(tempA, 3)); 	\
}

// CPU intrinsic popcount
// This variation appears to be the fastest in our application
// Most likely because a result vector is already in an anonymous register
// and does not need to be reloaded
//
// _popcnt64 performance on all architectures:
// Latency 3, Throughput 1
#define POPCOUNT2(A, B) { 							 	\
	__m256i tempA = _mm512_extracti64x4_epi64(B, 0);	\
	A += _popcnt64(_mm256_extract_epi64(tempA, 0));		\
	A += _popcnt64(_mm256_extract_epi64(tempA, 1)); 	\
	A += _popcnt64(_mm256_extract_epi64(tempA, 2)); 	\
	A += _popcnt64(_mm256_extract_epi64(tempA, 3)); 	\
	tempA = _mm512_extracti64x4_epi64(B, 1); 			\
	A += _popcnt64(_mm256_extract_epi64(tempA, 0));		\
	A += _popcnt64(_mm256_extract_epi64(tempA, 1)); 	\
	A += _popcnt64(_mm256_extract_epi64(tempA, 2)); 	\
	A += _popcnt64(_mm256_extract_epi64(tempA, 3)); 	\
}

#define FILTER_UNPHASED(A, B)			 _mm512_and_si512(_mm512_slli_epi64(_mm512_and_si512(_mm512_or_si512(_mm512_and_si512(A, maskUnphasedHigh),_mm512_and_si512(B, maskUnphasedLow)), maskUnphasedLow), 1), A)
#define FILTER_UNPHASED_PAIR(A, B, C, D) _mm512_or_si512(_mm512_srli_epi64(FILTER_UNPHASED(A, B), 1), FILTER_UNPHASED(C, D))
#define FILTER_UNPHASED_SPECIAL(A)		 _mm512_and_si512(_mm512_and_si512(_mm512_srli_epi64(A, 1), A), maskUnphasedLow)


#elif SIMD_VERSION == 5 // AVX2
#define VECTOR_TYPE	__m256i
const VECTOR_TYPE ONE_MASK = _mm256_set1_epi8(255); // 11111111b
const VECTOR_TYPE maskUnphasedHigh = _mm256_set1_epi8(UNPHASED_UPPER_MASK);	// 10101010b
const VECTOR_TYPE maskUnphasedLow  = _mm256_set1_epi8(UNPHASED_LOWER_MASK);	// 01010101b

#define PHASED_ALTALT(A,B)	_mm256_and_si256(A, B)
#define PHASED_REFREF(A,B)	_mm256_and_si256(_mm256_xor_si256(A, ONE_MASK), _mm256_xor_si256(B, ONE_MASK))
#define PHASED_ALTREF(A,B)	_mm256_and_si256(_mm256_xor_si256(A, B), B)
#define PHASED_REFALT(A,B)	_mm256_and_si256(_mm256_xor_si256(A, B), A)
#define PHASED_ALTALT_MASK(A,B,M)	_mm256_and_si256(PHASED_ALTALT(A, B), M)
#define PHASED_REFREF_MASK(A,B,M)	_mm256_and_si256(PHASED_REFREF(A, B), M)
#define PHASED_ALTREF_MASK(A,B,M)	_mm256_and_si256(PHASED_ALTREF(A, B), M)
#define PHASED_REFALT_MASK(A,B,M)	_mm256_and_si256(PHASED_REFALT(A, B), M)
#define MASK_MERGE(A,B)		_mm256_xor_si256(_mm256_or_si256(A, B), ONE_MASK)

// Software intrinsic popcount
#define POPCOUNT(A, B) {									\
	A += __builtin_popcountll(_mm256_extract_epi64(B, 0));	\
	A += __builtin_popcountll(_mm256_extract_epi64(B, 1));	\
	A += __builtin_popcountll(_mm256_extract_epi64(B, 2));	\
	A += __builtin_popcountll(_mm256_extract_epi64(B, 3));	\
}

// CPU intrinsic popcount
#define POPCOUNT2(A, B) { 						\
	A += _popcnt64(_mm256_extract_epi64(B, 0));	\
	A += _popcnt64(_mm256_extract_epi64(B, 1)); \
	A += _popcnt64(_mm256_extract_epi64(B, 2)); \
	A += _popcnt64(_mm256_extract_epi64(B, 3)); \
}

#define FILTER_UNPHASED(A, B)			 _mm256_and_si256(_mm256_slli_epi64(_mm256_and_si256(_mm256_or_si256(_mm256_and_si256(A, maskUnphasedHigh),_mm256_and_si256(B, maskUnphasedLow)), maskUnphasedLow), 1), A)
#define FILTER_UNPHASED_PAIR(A, B, C, D) _mm256_or_si256(_mm256_srli_epi64(FILTER_UNPHASED(A, B), 1), FILTER_UNPHASED(C, D))
#define FILTER_UNPHASED_SPECIAL(A)		 _mm256_and_si256(_mm256_and_si256(_mm256_srli_epi64(A, 1), A), maskUnphasedLow)

#elif SIMD_VERSION >= 2 // SSE2+
#define VECTOR_TYPE	__m128i
const VECTOR_TYPE ONE_MASK = _mm_set1_epi8(255); // 11111111b
const VECTOR_TYPE maskUnphasedHigh = _mm_set1_epi8(UNPHASED_UPPER_MASK);	// 10101010b
const VECTOR_TYPE maskUnphasedLow  = _mm_set1_epi8(UNPHASED_LOWER_MASK);	// 01010101b

#define PHASED_ALTALT(A,B)	_mm_and_si128(A, B)
#define PHASED_REFREF(A,B)	_mm_and_si128(_mm_xor_si128(A, ONE_MASK), _mm_xor_si128(B, ONE_MASK))
#define PHASED_ALTREF(A,B)	_mm_and_si128(_mm_xor_si128(A, B), B)
#define PHASED_REFALT(A,B)	_mm_and_si128(_mm_xor_si128(A, B), A)
#define PHASED_ALTALT_MASK(A,B,M)	_mm_and_si128(PHASED_ALTALT(A, B), M)
#define PHASED_REFREF_MASK(A,B,M)	_mm_and_si128(PHASED_REFREF(A, B), M)
#define PHASED_ALTREF_MASK(A,B,M)	_mm_and_si128(PHASED_ALTREF(A, B), M)
#define PHASED_REFALT_MASK(A,B,M)	_mm_and_si128(PHASED_REFALT(A, B), M)
#define MASK_MERGE(A,B)		_mm_xor_si128(_mm_or_si128(A, B), ONE_MASK)

#define POPCOUNT(A, B) {									\
	A += __builtin_popcountll(_mm_extract_epi64(B, 0));		\
	A += __builtin_popcountll(_mm_extract_epi64(B, 1));		\
}

#define POPCOUNT2(A, B) { 						\
	A += _popcnt64(_mm_extract_epi64(B, 0));	\
	A += _popcnt64(_mm_extract_epi64(B, 1)); 	\
}

#define FILTER_UNPHASED(A, B)			 _mm_and_si128(_mm_slli_epi64(_mm_and_si128(_mm_or_si128(_mm_and_si128(A, maskUnphasedHigh),_mm_and_si128(B, maskUnphasedLow)), maskUnphasedLow), 1), A)
#define FILTER_UNPHASED_PAIR(A, B, C, D) _mm_or_si128(_mm_srli_epi64(FILTER_UNPHASED(A, B), 1), FILTER_UNPHASED(C, D))
#define FILTER_UNPHASED_SPECIAL(A)		 _mm_and_si128(_mm_and_si128(_mm_srli_epi64(A, 1), A), maskUnphasedLow)
#endif
#endif // ENDIF SIMD_AVAILABLE == 1

template <class T>
class TomahawkCalculateSlave{
	//Basic typedefs
	typedef TomahawkCalculateSlave<T> self_type;
	typedef TomahawkReaderControllerManager<const T> manager_type;
	typedef TomahawkReaderController<const T> controller_type;
	typedef const TomahawkEntryMeta<const T> meta_type;
	typedef const Support::TomahawkRun<const T> run_type;
	typedef TotempoleEntry totempole_entry_type;
	typedef TomahawkCalculateSlaveOutputManager<T> output_manager_type;
	typedef IO::TomahawkCalculationWriterInterace writer_type;
	typedef Support::TomahawkOutputLD helper_type;

	typedef TomahawkReaderControllerPackedPair<> simd_pair;

	// Work orders
	typedef Tomahawk::LoadBalancerBlock order_type;
	typedef std::vector<order_type> work_order;

	// Function pointers
	typedef void (self_type::*phaseFunction)(const controller_type& block1, const controller_type block2);

public:
	TomahawkCalculateSlave(const manager_type& manager,
						   writer_type& writer,
						   Interface::ProgressBar& progress,
						   const TomahawkCalcParameters& parameters,
						   const work_order& orders);

	~TomahawkCalculateSlave();
	TomahawkCalculateSlave(const TomahawkCalculateSlave& other);
	TomahawkCalculateSlave(TomahawkCalculateSlave&& other) noexcept;
	TomahawkCalculateSlave& operator=(const TomahawkCalculateSlave& other);
	TomahawkCalculateSlave& operator=(TomahawkCalculateSlave&& other) noexcept;
	TomahawkCalculateSlave& operator+=(const TomahawkCalculateSlave& other);

	std::thread* Start(void){
		this->thread = std::thread(&self_type::Calculate, this);
		return(&this->thread);
	}

	inline const U64& getImpossible(void) const{ return this->impossible; }
	inline const U64& getPossible(void) const{ return this->possible; }
	inline const U64& getNoHets(void) const{ return this->no_uncertainty; }
	inline const U64& getInsufficientData(void) const{ return this->insufficent_alleles; }
	inline U64 getComparisons(void) const{ return(this->impossible + this->possible + this->insufficent_alleles); }

private:
	bool Calculate();

	bool DiagonalWorkOrder(const order_type& order);
	bool SquareWorkOrder(const order_type& order);
	bool CompareBlocks(controller_type& block1);
	bool CompareBlocks(controller_type& block1, controller_type block2);
	inline void CompareBlocksFunction(const controller_type& block1, const controller_type block2);
	inline void CompareBlocksFunctionForcedPhased(const controller_type& block1, const controller_type block2);
	inline void CompareBlocksFunctionForcedUnphased(const controller_type& block1, const controller_type block2);

	bool CalculateLDPhased(const controller_type& a, const controller_type& b);
	bool CalculateLDPhasedMath(void);
	bool CalculateLDPhasedVectorized(const controller_type& a, const controller_type& b);
	bool CalculateLDPhasedVectorizedNoMissing(const controller_type& a, const controller_type& b);
	bool CalculateLDUnphased(const controller_type& a, const controller_type& b);
	bool CalculateLDUnphasedVectorized(const controller_type& a, const controller_type& b);
	bool CalculateLDUnphasedVectorizedNoMissing(const controller_type& a, const controller_type& b);
	bool CalculateLDUnphasedMath(void);

	double EstimateChiSq(const double& target, const double& p, const double& q) const;
	bool ChooseF11Calculate(const double& target, const double& p, const double& q);

	void setFLAGs(const controller_type& a, const controller_type& b);

private:
	const TomahawkCalcParameters& parameters;

	// Counters
	U32 block_comparisons;
	U64 variant_comparisons;

	// Counters
	const U64 samples;
	U64 impossible;
	U64 possible;
	U64 no_uncertainty;
	U64 insufficent_alleles;
	//U64 false_positive;
	//U64 false_negative;

	helper_type helper;
	Support::TomahawkSlaveSIMDHelper<> helper_simd;
	Algorithm::FisherTest fisherController;
	const manager_type& manager;
	IO::BasicBuffer outstreamBuffer;

	// thread
	std::thread thread;

	// writer manager
	output_manager_type output_manager;

	// progress
	Interface::ProgressBar& progress;

	// function pointers
	phaseFunction phase_function_across;

	// thread work orders
	work_order orders;

	// SIMD parameters
	const U32 byte_width; // Number of bytes required per variant site
	const U32 byteAlignedEnd; // End byte position
	const U32 vectorCycles; // Number of SIMD cycles (genotypes/2/vector width)
	const U32 phased_unbalanced_adjustment; // Modulus remainder
	const U32 unphased_unbalanced_adjustment; // Modulus remainder
};

template <class T>
void TomahawkCalculateSlave<T>::setFLAGs(const controller_type& a, const controller_type& b){
	// If long range
	const meta_type& mA = a.currentMeta();
	const meta_type& mB = b.currentMeta();

	if(b.support->contigID == a.support->contigID)
		this->helper.setSameContig();

	if((mB.position >> 2) - (mA.position >> 2) > LONG_RANGE_THRESHOLD)
		this->helper.setLongRange();

	// Set FLAGs
	if(mA.MAF < LOW_MAF_THRESHOLD)
		this->helper.setLowMAFA();

	if(mB.MAF < LOW_MAF_THRESHOLD)
		this->helper.setLowMAFB();

	if(mA.HWE_P < LOW_HWE_THRESHOLD)
		this->helper.setFailedHWEA();

	if(mB.HWE_P < LOW_HWE_THRESHOLD)
		this->helper.setFailedHWEB();

	if(mA.missing || mB.missing)
		this->helper.setHasMissingValues();
}

template <class T>
TomahawkCalculateSlave<T>::TomahawkCalculateSlave(const manager_type& manager,
						   writer_type& writer,
						   Interface::ProgressBar& progress,
						   const TomahawkCalcParameters& parameters,
						   const work_order& orders) :
		parameters(parameters),
		block_comparisons(0),
		variant_comparisons(0),
		samples(manager.header.getSamples()),
		impossible(0),
		possible(0),
		no_uncertainty(0),
		insufficent_alleles(0),
		//false_positive(0),
		//false_negative(0),
		fisherController(1024),
		manager(manager),
		output_manager(writer, parameters.compression_type),
		progress(progress),
		phase_function_across(nullptr),
		orders(orders),
		byte_width(ceil((double)samples/4)),
		byteAlignedEnd(this->byte_width/(GENOTYPE_TRIP_COUNT/4)*(GENOTYPE_TRIP_COUNT/4)),
		vectorCycles(this->byteAlignedEnd*4/GENOTYPE_TRIP_COUNT),
		phased_unbalanced_adjustment((this->samples*2)%8),
		unphased_unbalanced_adjustment(this->samples%4)
	{
		if(this->parameters.force == TomahawkCalcParameters::force_method::none)
			this->phase_function_across = &self_type::CompareBlocksFunction;
		else if(this->parameters.force == TomahawkCalcParameters::force_method::phasedFunction)
			this->phase_function_across = &self_type::CompareBlocksFunctionForcedPhased;
		else
			this->phase_function_across = &self_type::CompareBlocksFunctionForcedUnphased;
	}

	template <class T>
	TomahawkCalculateSlave<T>::~TomahawkCalculateSlave(){ }

	/** Copy constructor */
	template <class T>
	TomahawkCalculateSlave<T>::TomahawkCalculateSlave (const TomahawkCalculateSlave& other) :
		parameters(other.parameters),
		block_comparisons(other.block_comparisons),
		variant_comparisons(other.variant_comparisons),
		samples(other.samples),
		impossible(other.impossible),
		possible(other.possible),
		no_uncertainty(other.no_uncertainty),
		insufficent_alleles(other.insufficent_alleles),
		//false_positive(other.false_positive),
		//false_negative(other.false_negative),
		fisherController(1024),
		manager(manager),
		output_manager(other.output_manager),
		progress(other.progress),
		phase_function_across(other.phase_function_across),
		orders(other.orders),
		byte_width(other.byte_width),
		byteAlignedEnd(other.byteAlignedEnd),
		vectorCycles(other.vectorCycles),
		phased_unbalanced_adjustment(other.phased_unbalanced_adjustment),
		unphased_unbalanced_adjustment(other.unphased_unbalanced_adjustment)
	{

	}

	/** Move constructor */
	template <class T>
	TomahawkCalculateSlave<T>::TomahawkCalculateSlave(TomahawkCalculateSlave&& other) noexcept :
		parameters(other.parameters),
		block_comparisons(other.block_comparisons),
		variant_comparisons(other.variant_comparisons),
		samples(other.samples),
		impossible(other.impossible),
		possible(other.possible),
		no_uncertainty(other.no_uncertainty),
		insufficent_alleles(other.insufficent_alleles),
		//false_positive(other.false_positive),
		//false_negative(other.false_negative),
		fisherController(1024),
		manager(other.manager),
		output_manager(other.output_manager),
		progress(other.progress),
		phase_function_across(other.phase_function_across),
		orders(other.orders),
		byte_width(other.byte_width),
		byteAlignedEnd(other.byteAlignedEnd),
		vectorCycles(other.vectorCycles),
		phased_unbalanced_adjustment(other.phased_unbalanced_adjustment),
		unphased_unbalanced_adjustment(other.unphased_unbalanced_adjustment)
	{
	   std::swap(this->thread, other.thread);
	}

	 /** Copy assignment operator */
	template <class T>
	TomahawkCalculateSlave<T>& TomahawkCalculateSlave<T>::operator=(const TomahawkCalculateSlave<T>& other){
		self_type tmp(other);         // re-use copy-constructor
		*this = std::move(tmp); // re-use move-assignment
		return *this;
	}

	/** Move assignment operator */
	template <class T>
	TomahawkCalculateSlave<T>& TomahawkCalculateSlave<T>::operator=(TomahawkCalculateSlave<T>&& other) noexcept{
		std::swap(this->thread, other.thread);
		return *this;
	}

	// Reduce function
	template <class T>
	TomahawkCalculateSlave<T>& TomahawkCalculateSlave<T>::operator+=(const TomahawkCalculateSlave<T>& other){
		this->block_comparisons += other.block_comparisons;
		this->variant_comparisons += other.variant_comparisons;
		this->impossible += other.impossible;
		this->possible += other.possible;
		this->no_uncertainty += other.no_uncertainty;
		this->insufficent_alleles += other.insufficent_alleles;
		//this->false_positive += other.false_positive;
		//this->false_negative += other.false_negative;

		return(*this);
	}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDUnphased(const controller_type& a, const controller_type& b){
	if(a.meta[a.metaPointer].MAF == 0 || b.meta[b.metaPointer].MAF == 0)
		return false;

	this->helper.resetUnphased();

	/////////////
	// Calculate
	/////////////
	T currentLengthA = a[0].runs;
	T currentLengthB = b[0].runs;

	BYTE currentMix = ((a[0].alleleA & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH*3)
					^ ((a[0].alleleB & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH*2)
					^ ((b[0].alleleA & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH*1)
					^ ((b[0].alleleB & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)));

	U32 pointerA = 0;
	U32 pointerB = 0;

#if __DEBUG_MODE__ == 4
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	while(true){
		//if(currentMix != 0 || currentMix != 1 || currentMix != 5 || currentMix != 16 || currentMix != 17 || currentMix != 21 || currentMix != 80 || currentMix != 81 || currentMix != 84)
		//std::cerr << (int)currentMix << std::endl;
		// If processed run length A > processed run length B
		if(currentLengthA > currentLengthB){
			currentLengthA -= currentLengthB;
			this->helper[currentMix] += currentLengthB;
			++pointerB;
			currentLengthB = b[pointerB].runs;
		}
		// If processed run length A == processed run length B
		else if(currentLengthA == currentLengthB){
			this->helper[currentMix] += currentLengthB;
			++pointerA;
			++pointerB;
			currentLengthA = a[pointerA].runs;
			currentLengthB = b[pointerB].runs;
		}
		// If processed run length A < processed run length B
		else {
			currentLengthB -= currentLengthA;
			this->helper[currentMix] += currentLengthA;
			++pointerA;
			currentLengthA = a[pointerA].runs;
		}

		// Exit condition
		if(pointerA == a.meta[a.metaPointer].runs || pointerB == b.meta[b.metaPointer].runs){
			//std::cerr << pointerA << '/' << a.meta[a.metaPointer].runs << '\t' << pointerB << '/' << b.meta[b.metaPointer].runs << std::endl;
			if(pointerA != a.meta[a.metaPointer].runs || pointerB != b.meta[b.metaPointer].runs){
				std::cerr << Tomahawk::Helpers::timestamp("FATAL") << "Failed to exit equally!\n" << pointerA << "/" << a.meta[a.metaPointer].runs << " and " << pointerB << "/" << b.meta[b.metaPointer].runs << std::endl;
				exit(1);
			}
			break;
		}

		currentMix  = ((a[pointerA].alleleA & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH*3)
					^ ((a[pointerA].alleleB & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH*2)
					^ ((b[pointerB].alleleA & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH*1)
					^ ((b[pointerB].alleleB & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)));
	}

#if __DEBUG_MODE__ == 4
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cerr << "T\t" << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\t'
				<< this->helper.alleleCounts[0] << '\t' << this->helper.alleleCounts[1] + this->helper.alleleCounts[4] << '\t' << this->helper.alleleCounts[5] << '\t'
				  << this->helper.alleleCounts[16] + this->helper.alleleCounts[64] << '\t' << this->helper.alleleCounts[17] + this->helper.alleleCounts[20] + this->helper.alleleCounts[65] + this->helper.alleleCounts[68] << '\t' << this->helper.alleleCounts[21] + this->helper.alleleCounts[69] << '\t'
				  << this->helper.alleleCounts[80] << '\t' << this->helper.alleleCounts[81]+this->helper.alleleCounts[84] << '\t' << this->helper.alleleCounts[85] << std::endl;
#elif __DEBUG_MODE__ == 6
	std::cerr << "U\t"
			<< this->helper.alleleCounts[0] << '\t' << this->helper.alleleCounts[1] + this->helper.alleleCounts[4] << '\t' << this->helper.alleleCounts[5] << '\t'
		  << this->helper.alleleCounts[16] + this->helper.alleleCounts[64] << '\t' << this->helper.alleleCounts[17] + this->helper.alleleCounts[20] + this->helper.alleleCounts[65] + this->helper.alleleCounts[68] << '\t' << this->helper.alleleCounts[21] + this->helper.alleleCounts[69] << '\t'
		  << this->helper.alleleCounts[80] << '\t' << this->helper.alleleCounts[81]+this->helper.alleleCounts[84] << '\t' << this->helper.alleleCounts[85] << std::endl;
#endif

	this->setFLAGs(a, b);

	return(this->CalculateLDUnphasedMath());
}

template <class T>
double TomahawkCalculateSlave<T>::EstimateChiSq(const double& target, const double& p, const double& q) const{
	const double f12 = p - target;
	const double f21 = q - target;
	const double f22 = 1 - (target + f12 + f21);
	const double e1111 = this->helper.totalAlleleCounts * pow(target,2);
	const double e1112 = 2 * this->helper.totalAlleleCounts * target * f12;
	const double e1122 = this->helper.totalAlleleCounts * pow(f12,2);
	const double e1211 = 2 * this->helper.totalAlleleCounts * target * f21;
	const double e1212 = 2 * this->helper.totalAlleleCounts * f12 * f21 + 2 * this->helper.totalAlleleCounts * target * f22;
	const double e1222 = 2 * this->helper.totalAlleleCounts * f12 * f22;
	const double e2211 = this->helper.totalAlleleCounts * pow(f21,2);
	const double e2212 = 2 * this->helper.totalAlleleCounts * f21 * f22;
	const double e2222 = this->helper.totalAlleleCounts * pow(f22,2);
	const double chisq1111 = e1111 > 0 ? pow(this->helper[0] - e1111, 2) / e1111 : 0,
				 chisq1112 = e1112 > 0 ? pow(this->helper[1] + this->helper[4] - e1112, 2) / e1112 : 0,
				 chisq1122 = e1122 > 0 ? pow(this->helper[5] - e1122, 2) / e1122 : 0,
				 chisq1211 = e1211 > 0 ? pow(this->helper[16] + this->helper[64] - e1211, 2) / e1211 : 0,
				 chisq1212 = e1212 > 0 ? pow(this->helper[17] + this->helper[20] + this->helper[65] + this->helper[68] - e1212, 2) / e1212 : 0,
				 chisq1222 = e1222 > 0 ? pow(this->helper[21] + this->helper[69] - e1222, 2) / e1222 : 0,
				 chisq2211 = e2211 > 0 ? pow(this->helper[80] - e2211, 2) / e2211 : 0,
				 chisq2212 = e2212 > 0 ? pow(this->helper[81] + this->helper[84] - e2212, 2) / e2212 : 0,
				 chisq2222 = e2222 > 0 ? pow(this->helper[85] - e2222, 2) / e2222 : 0;

	return(chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222);
}

template <class T>
bool TomahawkCalculateSlave<T>::ChooseF11Calculate(const double& target, const double& p, const double& q){
	this->helper.haplotypeCounts[0] = target;
	this->helper.haplotypeCounts[1] = p - this->helper.haplotypeCounts[0];
	this->helper.haplotypeCounts[2] = q - this->helper.haplotypeCounts[0];
	this->helper.haplotypeCounts[3] = 1 - (this->helper.haplotypeCounts[0] + this->helper.haplotypeCounts[1] + this->helper.haplotypeCounts[2]);

	this->helper.D = this->helper.haplotypeCounts[0] * this->helper.haplotypeCounts[3] - this->helper.haplotypeCounts[1] * this->helper.haplotypeCounts[2];
	this->helper.R2 = this->helper.D*this->helper.D / (p * (1 - p) * q * (1 - q));

	if(this->helper.R2 >= this->parameters.R2_min && this->helper.R2 <= this->parameters.R2_max){
		if(this->helper.D >= 0){
			this->helper.Dmax = p * (1.0-q) < q * (1.0-p)
					? p * (1.0-q)
					: q * (1.0-p);
		} else {
			this->helper.Dmax = p*q < (1-p)*(1-q)
					? p*q
					: (1-p)*(1-q);
		}
		this->helper.Dprime = this->helper.D / this->helper.Dmax;

		// Calculate p, q here
		// estimate p,q
		const double e1111 = this->helper.totalAlleleCounts * pow(this->helper.haplotypeCounts[0],2);
		const double e1112 = 2 * this->helper.totalAlleleCounts * this->helper.haplotypeCounts[0] * this->helper.haplotypeCounts[1];
		const double e1122 = this->helper.totalAlleleCounts * pow(this->helper.haplotypeCounts[1],2);
		const double e1211 = 2 * this->helper.totalAlleleCounts * this->helper.haplotypeCounts[0] * this->helper.haplotypeCounts[2];
		const double e1212 = 2 * this->helper.totalAlleleCounts * this->helper.haplotypeCounts[1] * this->helper.haplotypeCounts[2] + 2 * this->helper.totalAlleleCounts * this->helper.haplotypeCounts[0] * this->helper.haplotypeCounts[3];
		const double e1222 = 2 * this->helper.totalAlleleCounts * this->helper.haplotypeCounts[1] * this->helper.haplotypeCounts[3];
		const double e2211 = this->helper.totalAlleleCounts * pow(this->helper.haplotypeCounts[2],2);
		const double e2212 = 2 * this->helper.totalAlleleCounts * this->helper.haplotypeCounts[2] * this->helper.haplotypeCounts[3];
		const double e2222 = this->helper.totalAlleleCounts * pow(this->helper.haplotypeCounts[3],2);
		// const double total = e1111 + e1112 + e1122 + e1211 + e1212 + e1222 + e2211 + e2212 + e2222;

		const double p1 = 2*e1111 + e1112 + e1211 + e1212/4;
		const double p2 = e1112 + 2*e1122 + e1212/4 + e1222;
		const double q1 = e1211 + e1212/4 + 2*e2211 + e2212;
		const double q2 = e1212/4 + e1222 + e2212 + 2*e2222;

		if(this->helper.countAlternatives() < this->parameters.minimum_alleles)
			return false;

		// Calculate P: Fisher's exact test
		this->helper.chiSqFisher = this->fisherController.chiSquaredTest(p1,p2,q1,q2);

		if(this->helper.chiSqFisher > CHI_SQ_MAX_CV){ // Rough lower limit for CV that is possible to calculate
			const U64 a = ceil(p1);
			const U64 b = ceil(p2);
			const U64 c = ceil(q1);
			const U64 d = ceil(q2);

			// Toggle Fisher's exact test FLAG bit
			this->helper.setFisherTest();

			if(this->helper.D < 0)
				this->helper.P = this->fisherController.fisherTestLess(a,b,c,d);
			else
				this->helper.P = this->fisherController.fisherTestGreater(a,b,c,d);
		} else
			this->helper.P = this->fisherController.chisqr(1, this->helper.chiSqFisher);

		// Fisher's exact test P value filter
		if(this->helper.P > this->parameters.P_threshold)
			return false;

		// For writing output
		this->helper[0] = p1;
		this->helper[1] = p2;
		this->helper[4] = q1;
		this->helper[5] = q2;

		return true;
	}

	return false;
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDUnphasedMath(void){
	// Total amount of non-missing alleles
	this->helper.totalAlleleCounts  = this->helper[0]  + this->helper[1]  + this->helper[4]  + this->helper[5]
									+ this->helper[16] + this->helper[17] + this->helper[20] + this->helper[21]
									+ this->helper[64] + this->helper[65] + this->helper[68] + this->helper[69]
									+ this->helper[80] + this->helper[81] + this->helper[84] + this->helper[85];

	// All values are missing or too few
	if(this->helper.totalAlleleCounts < MINIMUM_ALLOWED_ALLELES){
		++this->insufficent_alleles;
		return false;
	}

	// How many hets-hets is there? 0/1 0/1 or 1/0 1/0 or equivalent
	const float number_of_hets = this->helper[17] + this->helper[20] + this->helper[65] + this->helper[68];

	// If het-hets are 0 then do normal calculations
	// There is no phase uncertainty
	// Use phased math
	if(number_of_hets == 0){
		this->helper[0] = 2*this->helper[0] + this->helper[1] + this->helper[4] + this->helper[16] + this->helper[17] + this->helper[64] + this->helper[68];
		this->helper[1] = this->helper[1] + this->helper[4] + 2*this->helper[5] + this->helper[20] + this->helper[21] + this->helper[65] + this->helper[69];
		this->helper[4] = this->helper[16] + this->helper[20] + this->helper[64] + this->helper[65] + 2*this->helper[80] + this->helper[81] + this->helper[84];
		this->helper[5] = this->helper[17] + this->helper[21] + this->helper[68] + this->helper[69] + this->helper[81] + this->helper[84] + 2*this->helper[85];

		++this->no_uncertainty;

		// Use standard phased LD math
		return(this->CalculateLDPhasedMath());
	}

	const double p = ((this->helper[0] + this->helper[1] + this->helper[4] + this->helper[5])*2.0
				   + (this->helper[16] + this->helper[17] + this->helper[20] + this->helper[21] + this->helper[64] + this->helper[65] + this->helper[68] + this->helper[69]))
				   / (2.0 * this->helper.totalAlleleCounts);
	const double q = ((this->helper[0] + this->helper[16] + this->helper[64] + this->helper[80])*2.0
				   + (this->helper[1] + this->helper[4] + this->helper[17] + this->helper[20] + this->helper[65] + this->helper[68] + this->helper[81] + this->helper[84]))
				   / (2.0 * this->helper.totalAlleleCounts);
	const double n11 = (2.0* this->helper[0] + this->helper[1] + this->helper[4] + this->helper[16] + this->helper[64]);

	// Not used for anything
	//const double n12 = (2.0*this->helper[5]  + this->helper[1]  + this->helper[4]  + this->helper[21] + this->helper[69]);
	//const double n21 = (2.0*this->helper[80] + this->helper[81] + this->helper[84] + this->helper[16] + this->helper[64]);
	//const double n22 = (2.0*this->helper[85] + this->helper[81] + this->helper[84] + this->helper[21] + this->helper[85]);

	/////////////////////////
	// Cubic function: a3x^3 + a2x^2 + a1x + d = 0 <==> ax^3 + bx^2 + cx + d = 0
	// Cubic constants
	/////////////////////////
	const double G   = 1.0 - 2.0*p - 2.0*q;
	const double dee = -n11*p*q;
	const double c   = -n11*G - number_of_hets*(1.0 - p - q) + 2.0*this->helper.totalAlleleCounts*p*q;
	const double b   = 2.0*this->helper.totalAlleleCounts*G - 2.0*n11 - number_of_hets;
	const double a   = 4.0 * this->helper.totalAlleleCounts;

	// Bounds for biological relevance
	const double minhap = n11 / (2.0 * this->helper.totalAlleleCounts);
	const double maxhap = (n11 + number_of_hets) / (2.0 * this->helper.totalAlleleCounts);

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
		const double theta 		= acos(-yN / sqrt(h2)) / 3.0;
		const double constant 	= 2.0 * sqrt(d2);
		const double alpha 		= xN + constant * cos(theta);
		const double beta 		= xN + constant * cos(2.0*M_PI/3.0 + theta);
		const double gamma 		= xN + constant * cos(4.0*M_PI/3.0 + theta);

		BYTE biologically_possible = 0;
		this->helper.chiSqModel = std::numeric_limits<float>::max();
		const double* chosen = &alpha;
		if(alpha >= minhap - Constants::ALLOWED_ROUNDING_ERROR && alpha <= maxhap + Constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			this->helper.chiSqModel = this->EstimateChiSq(alpha, p, q);
			chosen = &alpha;
		}

		if(beta >= minhap - Constants::ALLOWED_ROUNDING_ERROR && beta <= maxhap + Constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			if(this->EstimateChiSq(beta, p, q) < this->helper.chiSqModel){
				chosen = &beta;
				this->helper.chiSqModel = this->EstimateChiSq(beta, p, q);
			}
		}

		if(gamma >= minhap - Constants::ALLOWED_ROUNDING_ERROR && gamma <= maxhap + Constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			if(this->EstimateChiSq(gamma, p, q) < this->helper.chiSqModel){
				chosen = &gamma;
				this->helper.chiSqModel = this->EstimateChiSq(gamma, p, q); // implicit
			}
		}

		if(biologically_possible == 0){
			++this->impossible;
			return false;
		}

		if(biologically_possible > 1)
			this->helper.setMultipleRoots();

		++this->possible;
		return(this->ChooseF11Calculate(*chosen, p, q));

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
		if(!(alpha >= minhap - Constants::ALLOWED_ROUNDING_ERROR && alpha <= maxhap + Constants::ALLOWED_ROUNDING_ERROR)){
			++this->impossible;
			return false;
		}

		this->helper.chiSqModel = this->EstimateChiSq(alpha, p, q);
		++this->possible;

		return(this->ChooseF11Calculate(alpha, p, q));

	} else { // Yn2 == h2
		const double delta = pow((yN/2.0*a),(1.0/3.0));
		const double alpha = xN + delta; // alpha = beta in this case
		const double gamma = xN - 2.0*delta;

		if(std::isnan(alpha) || std::isnan(gamma)){
			++this->impossible;
			return false;
		}

		BYTE biologically_possible = 0;
		this->helper.chiSqModel = std::numeric_limits<float>::max();
		const double* chosen = &alpha;
		if(alpha >= minhap - Constants::ALLOWED_ROUNDING_ERROR && alpha <= maxhap + Constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			this->helper.chiSqModel = this->EstimateChiSq(alpha, p, q);
			chosen = &alpha;
		}

		if(gamma >= minhap - Constants::ALLOWED_ROUNDING_ERROR && gamma <= maxhap + Constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			if(this->EstimateChiSq(gamma, p, q) < this->helper.chiSqModel){
				chosen = &gamma;
				this->helper.chiSqModel = this->EstimateChiSq(gamma, p, q); // implicit
			}
		}

		if(biologically_possible == 0){
			++this->impossible;
			return false;
		}

		++this->possible;
		return(this->ChooseF11Calculate(*chosen, p, q));
	}
	return(false);
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDUnphasedVectorizedNoMissing(const controller_type& a, const controller_type& b){
	this->helper.resetUnphased();

	this->helper_simd.counters[0] = 0;
	this->helper_simd.counters[1] = 0;
	this->helper_simd.counters[2] = 0;
	this->helper_simd.counters[3] = 0;
	this->helper_simd.counters[4] = 0;
	this->helper_simd.counters[5] = 0;
	this->helper_simd.counters[6] = 0;
	this->helper_simd.counters[7] = 0;
	this->helper_simd.counters[8] = 0;
	this->helper_simd.counters[9] = 0;

	const simd_pair& datA = a.packed->getData(a.metaPointer);
	const simd_pair& datB = b.packed->getData(b.metaPointer);
	const BYTE* const arrayA = datA.data;
	const BYTE* const arrayB = datB.data;

#if SIMD_AVAILABLE == 1
	const U32 frontSmallest = datA.frontZero < datB.frontZero ? datA.frontZero : datB.frontZero;
	const U32 tailSmallest = datA.tailZero < datB.tailZero ? datA.tailZero : datB.tailZero;
	U32 i = frontSmallest;
	const U32 frontBonus = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus = (datA.tailZero != tailSmallest ? datA.tailZero : datB.tailZero);

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	VECTOR_TYPE altalt, refref, altref, refalt;
	VECTOR_TYPE __intermediate;

// Debug timings
#if __DEBUG_MODE__ == 4
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
	__intermediate = FILTER_UNPHASED_SPECIAL(refref); 						\
	POPCOUNT2(this->helper_simd.counters[0], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);	\
	POPCOUNT2(this->helper_simd.counters[1], __intermediate);				\
	__intermediate = FILTER_UNPHASED(altref, altref); 						\
	POPCOUNT2(this->helper_simd.counters[2], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);	\
	POPCOUNT2(this->helper_simd.counters[3], __intermediate);				\
	__intermediate = FILTER_UNPHASED(refalt,refalt);						\
	POPCOUNT2(this->helper_simd.counters[6], __intermediate);				\
}

#define ITER_LONG {															\
	ITER_SHORT																\
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref); 	\
	POPCOUNT2(this->helper_simd.counters[5], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt); 	\
	POPCOUNT2(this->helper_simd.counters[7], __intermediate);				\
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt); 						\
	POPCOUNT2(this->helper_simd.counters[8], __intermediate);				\
}

	for( ; i < frontBonus; ) 					  	ITER_SHORT
	for( ; i < this->vectorCycles - tailBonus; )  	ITER_LONG
	for( ; i < this->vectorCycles - tailSmallest; ) ITER_SHORT

#undef ITER_LONG
#undef ITER_SHORT
#undef ITER_BASE

	U32 k = this->byteAlignedEnd;
#else
	U32 k = 0;
#endif

	BYTE b_altalt, b_refref, b_refalt, b_altref;
#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k < this->byte_width; ++k){
		b_altalt  = (arrayA[k] & arrayB[k]);
		b_refref  = ((~arrayA[k]) & (~arrayB[k]));
		b_altref  = ((arrayA[k] ^ arrayB[k]) & arrayA[k]);
		b_refalt  = ((arrayA[k] ^ arrayB[k]) & arrayB[k]);

		this->helper_simd.counters[0] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_SPECIAL(b_refref));
		this->helper_simd.counters[1] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_refalt, b_refalt, b_refref));
		this->helper_simd.counters[2] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE(b_refalt, b_refalt));
		this->helper_simd.counters[3] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_altref, b_altref, b_refref));
		//this->helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_altalt, b_altalt, b_refref));
		//this->helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refalt, b_altref, b_altref, b_refalt));

		this->helper_simd.counters[5] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refalt, b_altalt, b_altalt, b_refalt));
		this->helper_simd.counters[6] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE(b_altref, b_altref));
		this->helper_simd.counters[7] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_altref, b_altalt, b_altalt, b_altref));
		this->helper_simd.counters[8] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_SPECIAL(b_altalt));
	}

	this->helper[0]  = this->helper_simd.counters[0] - this->unphased_unbalanced_adjustment;
	this->helper[0] += (frontSmallest + tailSmallest) * GENOTYPE_TRIP_COUNT;
	this->helper[1]  = this->helper_simd.counters[1];
	this->helper[5]  = this->helper_simd.counters[2];
	this->helper[16] = this->helper_simd.counters[3];
	this->helper[21] = this->helper_simd.counters[5];
	this->helper[80] = this->helper_simd.counters[6];
	this->helper[81] = this->helper_simd.counters[7];
	this->helper[85] = this->helper_simd.counters[8];
	//this->helper[17] = this->helper_simd.counters[4];
	this->helper[17] = this->samples - (this->helper[0] +  this->helper[1] + this->helper[5] + this->helper[16] + this->helper[21] + this->helper[80] + this->helper[81] + this->helper[85]);

#if __DEBUG_MODE__ == 4
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cerr << "V\t" << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\t'
			  << this->helper.alleleCounts[0]  << '\t' << this->helper.alleleCounts[1]  << '\t' << this->helper.alleleCounts[5]  << '\t'
			  << this->helper.alleleCounts[16] << '\t' << this->helper.alleleCounts[17] << '\t' << this->helper.alleleCounts[21] << '\t'
			  << this->helper.alleleCounts[80] << '\t' << this->helper.alleleCounts[81] << '\t' << this->helper.alleleCounts[85] << std::endl;
#elif __DEBUG_MODE__ == 6
	std::cerr << "UVC\t"
			  << this->helper.alleleCounts[0]  << '\t' << this->helper.alleleCounts[1]  << '\t' << this->helper.alleleCounts[5]  << '\t'
			  << this->helper.alleleCounts[16] << '\t' << this->helper.alleleCounts[17] << '\t' << this->helper.alleleCounts[21] << '\t'
			  << this->helper.alleleCounts[80] << '\t' << this->helper.alleleCounts[81] << '\t' << this->helper.alleleCounts[85] << std::endl;
#endif

	this->setFLAGs(a, b);
	return(this->CalculateLDUnphasedMath());
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDUnphasedVectorized(const controller_type& a, const controller_type& b){
	#if __DEBUG_MODE__ < 6
	if(a.currentMeta().missing == 0 && b.currentMeta().missing == 0)
		return(this->CalculateLDUnphasedVectorizedNoMissing(a, b));
	#endif

	this->helper.resetUnphased();

	this->helper_simd.counters[0] = 0;
	this->helper_simd.counters[1] = 0;
	this->helper_simd.counters[2] = 0;
	this->helper_simd.counters[3] = 0;
	this->helper_simd.counters[4] = 0;
	this->helper_simd.counters[5] = 0;
	this->helper_simd.counters[6] = 0;
	this->helper_simd.counters[7] = 0;
	this->helper_simd.counters[8] = 0;
	this->helper_simd.counters[9] = 0;

	const simd_pair& datA = a.packed->getData(a.metaPointer);
	const simd_pair& datB = b.packed->getData(b.metaPointer);
	const BYTE* const arrayA = datA.data;
	const BYTE* const arrayB = datB.data;
	const BYTE* const arrayA_mask = datA.mask;
	const BYTE* const arrayB_mask = datB.mask;

#if SIMD_AVAILABLE == 1
	const U32 frontSmallest = datA.frontZero < datB.frontZero ? datA.frontZero : datB.frontZero;
	const U32 tailSmallest = datA.tailZero < datB.tailZero ? datA.tailZero : datB.tailZero;
	U32 i = frontSmallest;
	const U32 frontBonus = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus = (datA.tailZero != tailSmallest ? datA.tailZero : datB.tailZero);

	//std::cerr << frontSmallest << '\t' << tailSmallest << std::endl;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	const VECTOR_TYPE* const vectorA_mask = (const VECTOR_TYPE* const)arrayA_mask;
	const VECTOR_TYPE* const vectorB_mask = (const VECTOR_TYPE* const)arrayB_mask;
	VECTOR_TYPE altalt, refref, altref, refalt;
	//VECTOR_TYPE t00, t05, t50, t55, combTop, combLeft, combRight, combBottom;
	VECTOR_TYPE __intermediate, mask;

// Debug timings
#if __DEBUG_MODE__ == 4
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_BASE {												\
	mask	= MASK_MERGE(vectorA_mask[i], vectorB_mask[i]); 	\
	refref  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], mask);	\
	refalt  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], mask);	\
	altref  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], mask);	\
	altalt  = PHASED_ALTALT_MASK(vectorA[i], vectorB[i], mask);	\
	i += 1;														\
}

#define ITER_SHORT {														\
	ITER_BASE																\
	__intermediate = FILTER_UNPHASED_SPECIAL(refref); 						\
	POPCOUNT2(this->helper_simd.counters[0], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);	\
	POPCOUNT2(this->helper_simd.counters[1], __intermediate);				\
	__intermediate = FILTER_UNPHASED(altref, altref); 						\
	POPCOUNT2(this->helper_simd.counters[2], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);	\
	POPCOUNT2(this->helper_simd.counters[3], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref, altalt, altalt, refref);	\
	POPCOUNT2(this->helper_simd.counters[4], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altref, altref, refalt);	\
	POPCOUNT2(this->helper_simd.counters[4], __intermediate);				\
	__intermediate = FILTER_UNPHASED(refalt,refalt);						\
	POPCOUNT2(this->helper_simd.counters[6], __intermediate);				\
}

#define ITER_LONG {															\
	ITER_SHORT																\
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref); 	\
	POPCOUNT2(this->helper_simd.counters[5], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt); 	\
	POPCOUNT2(this->helper_simd.counters[7], __intermediate);				\
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt); 						\
	POPCOUNT2(this->helper_simd.counters[8], __intermediate);				\
}

	for( ; i < frontBonus; ) 					  	ITER_SHORT
	for( ; i < this->vectorCycles - tailBonus; )  	ITER_LONG
	for( ; i < this->vectorCycles - tailSmallest; ) ITER_SHORT

#undef ITER_LONG
#undef ITER_SHORT
#undef ITER_BASE
	U32 k = this->byteAlignedEnd;
#else
	U32 k = 0;
#endif

	BYTE b_altalt, b_refref, b_refalt, b_altref;
#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k < this->byte_width; ++k){
		b_altalt  = (arrayA[k] & arrayB[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_refref  = ((~arrayA[k]) & (~arrayB[k])) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_altref  = ((arrayA[k] ^ arrayB[k]) & arrayA[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_refalt  = ((arrayA[k] ^ arrayB[k]) & arrayB[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);

		this->helper_simd.counters[0] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_SPECIAL(b_refref));
		this->helper_simd.counters[1] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_refalt, b_refalt, b_refref));
		this->helper_simd.counters[2] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE(b_refalt, b_refalt));
		this->helper_simd.counters[3] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_altref, b_altref, b_refref));
		this->helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_altalt, b_altalt, b_refref));
		this->helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refalt, b_altref, b_altref, b_refalt));
		this->helper_simd.counters[5] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refalt, b_altalt, b_altalt, b_refalt));
		this->helper_simd.counters[6] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE(b_altref, b_altref));
		this->helper_simd.counters[7] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_altref, b_altalt, b_altalt, b_altref));
		this->helper_simd.counters[8] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_SPECIAL(b_altalt));
	}

	this->helper[0]  = this->helper_simd.counters[0] - this->unphased_unbalanced_adjustment;
	this->helper[0] += (frontSmallest + tailSmallest) * GENOTYPE_TRIP_COUNT;
	this->helper[1]  = this->helper_simd.counters[1];
	this->helper[5]  = this->helper_simd.counters[2];
	this->helper[16] = this->helper_simd.counters[3];
	this->helper[17] = this->helper_simd.counters[4];
	this->helper[21] = this->helper_simd.counters[5];
	this->helper[80] = this->helper_simd.counters[6];
	this->helper[81] = this->helper_simd.counters[7];
	this->helper[85] = this->helper_simd.counters[8];

#if __DEBUG_MODE__ == 4
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cerr << "V\t" << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\t'
			  << this->helper.alleleCounts[0]  << '\t' << this->helper.alleleCounts[1]  << '\t' << this->helper.alleleCounts[5]  << '\t'
			  << this->helper.alleleCounts[16] << '\t' << this->helper.alleleCounts[17] << '\t' << this->helper.alleleCounts[21] << '\t'
			  << this->helper.alleleCounts[80] << '\t' << this->helper.alleleCounts[81] << '\t' << this->helper.alleleCounts[85] << std::endl;
#elif __DEBUG_MODE__ == 6
	std::cerr << "UVM\t"
			  << this->helper.alleleCounts[0]  << '\t' << this->helper.alleleCounts[1]  << '\t' << this->helper.alleleCounts[5]  << '\t'
			  << this->helper.alleleCounts[16] << '\t' << this->helper.alleleCounts[17] << '\t' << this->helper.alleleCounts[21] << '\t'
			  << this->helper.alleleCounts[80] << '\t' << this->helper.alleleCounts[81] << '\t' << this->helper.alleleCounts[85] << std::endl;
#endif

	/*
	if(this->helper[0]+this->helper[1]+this->helper[5] + this->helper[16]+this->helper[17]+this->helper[21] + this->helper[80]+this->helper[81]+this->helper[85] != this->samples){
		std::cerr << a.metaPointer*a.packed->width << '/' << a.packed->total_size << '\t' << b.metaPointer*b.packed->width << '/' << b.packed->total_size << "\tExpected: " << this->samples << '/' << this->helper[0]+this->helper[1]+this->helper[5] + this->helper[16]+this->helper[17]+this->helper[21] + this->helper[80]+this->helper[81]+this->helper[85] << std::endl;
		std::cerr << this->manager[0].packed->total_size << '\t' << this->manager[1].packed->total_size << '\t' << this->manager[2].packed->total_size << '\t' << this->manager[3].packed->total_size << std::endl;
		exit(1);
	}
	*/

	this->setFLAGs(a, b);
	return(this->CalculateLDUnphasedMath());
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDPhasedVectorized(const controller_type& a, const controller_type& b){
	#if __DEBUG_MODE__ < 6
	if(a.currentMeta().missing == 0 && b.currentMeta().missing == 0)
		return(this->CalculateLDPhasedVectorizedNoMissing(a, b));
	#endif

	this->helper.resetPhased();
	this->helper_simd.counters[0] = 0;
	this->helper_simd.counters[1] = 0;
	this->helper_simd.counters[2] = 0;
	this->helper_simd.counters[3] = 0;

	const simd_pair& datA = a.packed->getData(a.metaPointer);
	const simd_pair& datB = b.packed->getData(b.metaPointer);
	const BYTE* const arrayA = datA.data;
	const BYTE* const arrayB = datB.data;
	const BYTE* const arrayA_mask = datA.mask;
	const BYTE* const arrayB_mask = datB.mask;

#if SIMD_AVAILABLE == 1
	const U32 frontSmallest = datA.frontZero < datB.frontZero ? datA.frontZero : datB.frontZero;
	const U32 tailSmallest = datA.tailZero < datB.tailZero ? datA.tailZero : datB.tailZero;
	U32 i = frontSmallest;
	const U32 frontBonus = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus = (datA.tailZero != tailSmallest ? datA.tailZero : datB.tailZero);

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	const VECTOR_TYPE* const vectorA_mask = (const VECTOR_TYPE* const)arrayA_mask;
	const VECTOR_TYPE* const vectorB_mask = (const VECTOR_TYPE* const)arrayB_mask;
	VECTOR_TYPE __intermediate, masks;

// Debug timings
#if __DEBUG_MODE__ == 4
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_SHORT {														\
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);					\
	__intermediate  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT2(this->helper_simd.counters[0], __intermediate);				\
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT2(this->helper_simd.counters[1], __intermediate);				\
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT2(this->helper_simd.counters[2], __intermediate);				\
	i += 1;																	\
}

#define ITER {																\
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);					\
	__intermediate  = PHASED_ALTALT_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT2(this->helper_simd.counters[3], __intermediate);				\
	__intermediate  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT2(this->helper_simd.counters[0], __intermediate);				\
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT2(this->helper_simd.counters[1], __intermediate);				\
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT2(this->helper_simd.counters[2], __intermediate);				\
	i += 1;																	\
}

	for( ; i < frontBonus; ) 					  	ITER_SHORT // Not possible to be ALT-ALT
	for( ; i < this->vectorCycles - tailBonus; )  	ITER
	for( ; i < this->vectorCycles - tailSmallest; ) ITER_SHORT // Not possible to be ALT-ALT

#undef ITER
#undef ITER_SHORT
	U32 k = this->byteAlignedEnd;
#else
	U32 k = 0;
#endif

	//std::cerr << datA.frontZero << '\t' << datB.frontZero << '\t' << datA.tailZero << '\t' << datB.tailZero << "\t->\t" << frontSmallest << '\t' << tailSmallest << std::endl;
	//std::cerr << frontSmallest << "->" << frontBonus << "->" << this->vectorCycles - tailBonus << "->" << this->vectorCycles - tailSmallest << "(" << (this->vectorCycles - tailSmallest)*(256/8) << ")/" << (this->byteAlignedEnd/(256/8)) << " skip " << this->byteAlignedEnd << "->" << this->byte_width << std::endl;

	BYTE mask;
#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k+8 < this->byte_width; k += 8){
#ifdef __INTEL_COMPILER
		#pragma vector aligned
#endif
		for(U32 l = 0; l < 8; ++l){
			mask = ~(arrayA_mask[k+l] | arrayB_mask[k+l]);
			this->helper_simd.scalarA[l] = (arrayA[k+l] & arrayB[k+l]) & mask;
			this->helper_simd.scalarB[l] = ((~arrayA[k+l]) & (~arrayB[k+l])) & mask;
			this->helper_simd.scalarC[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayA[k+l]) & mask;
			this->helper_simd.scalarD[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayB[k+l]) & mask;
		}
		this->helper_simd.counters[0] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(this->helper_simd.scalarB));
		this->helper_simd.counters[2] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(this->helper_simd.scalarC));
		this->helper_simd.counters[1] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(this->helper_simd.scalarD));
		this->helper_simd.counters[3] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(this->helper_simd.scalarA));
	}

#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k < this->byte_width; ++k){
		mask = ~(arrayA_mask[k] | arrayB_mask[k]);
		this->helper_simd.counters[0] += POPCOUNT_ITER(((~arrayA[k]) & (~arrayB[k])) & mask);
		this->helper_simd.counters[2] += POPCOUNT_ITER(((arrayA[k] ^ arrayB[k]) & arrayA[k]) & mask);
		this->helper_simd.counters[1] += POPCOUNT_ITER(((arrayA[k] ^ arrayB[k]) & arrayB[k]) & mask);
		this->helper_simd.counters[3] += POPCOUNT_ITER((arrayA[k] & arrayB[k]) & mask);
	}

	this->helper[1] = this->helper_simd.counters[1];
	this->helper[4] = this->helper_simd.counters[2];
	this->helper[5] = this->helper_simd.counters[3];
	//const float test = this->helper[0]; // Todo: fix to correct for missing values
	//this->helper[0] = this->samples*2 - (this->helper[1] + this->helper[4] + this->helper[5] + test + this->phased_unbalanced_adjustment); // Todo: have to correct this given missingness
	this->helper[0] = (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 + this->helper_simd.counters[0] - this->phased_unbalanced_adjustment;


#if __DEBUG_MODE__ == 4
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cerr << "V\t" << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << '\t' << this->helper[0] << '\t' << this->helper[1] << '\t' << this->helper[4] << '\t' << this->helper[5] << "\t" << this->helper[0]+this->helper[1]+this->helper[4]+this->helper[5] << "\t" << ticks_per_iter.count() << std::endl;
#elif __DEBUG_MODE__ == 6
	std::cerr << "PVM\t"
			  << this->helper[0] << '\t'
			  << this->helper[1] << '\t'
			  << this->helper[4] << '\t'
			  << this->helper[5] << std::endl;
#endif

	/*
	if(this->helper[0]+this->helper[1]+this->helper[4]+this->helper[5] != this->samples*2){
		std::cerr << "V\t" << this->helper[0] << '\t' << this->helper[1] << '\t' << this->helper[4] << '\t' << this->helper[5] << "\tsum: " << this->helper[0]+this->helper[1]+this->helper[4]+this->helper[5] << "\ttest: " << testSum << std::endl;
		std::cerr << a.metaPointer*a.packed->width << '/' << a.packed->total_size << '\t' << b.metaPointer*b.packed->width << '/' << b.packed->total_size << std::endl;
		std::cerr << this->manager[0].packed->total_size << '\t' << this->manager[1].packed->total_size << '\t' << this->manager[2].packed->total_size << '\t' << this->manager[3].packed->total_size << std::endl;
		exit(1);
	}
	*/

	this->setFLAGs(a, b);
	return(this->CalculateLDPhasedMath());
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDPhasedVectorizedNoMissing(const controller_type& a, const controller_type& b){
	this->helper.resetPhased();
	this->helper_simd.counters[0] = 0;
	this->helper_simd.counters[1] = 0;
	this->helper_simd.counters[2] = 0;
	this->helper_simd.counters[3] = 0;

	const simd_pair& datA = a.packed->getData(a.metaPointer);
	const simd_pair& datB = b.packed->getData(b.metaPointer);
	const BYTE* const arrayA = datA.data;
	const BYTE* const arrayB = datB.data;

#if SIMD_AVAILABLE == 1
	const U32 frontSmallest = datA.frontZero < datB.frontZero ? datA.frontZero : datB.frontZero;
	const U32 tailSmallest = datA.tailZero < datB.tailZero ? datA.tailZero : datB.tailZero;
	U32 i = frontSmallest;
	const U32 frontBonus = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus = (datA.tailZero != tailSmallest ? datA.tailZero : datB.tailZero);

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	//VECTOR_TYPE altalt, altref, refalt;
	VECTOR_TYPE __intermediate;

// Debug timings
#if __DEBUG_MODE__ == 4
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_SHORT {											\
	__intermediate  = PHASED_REFALT(vectorA[i], vectorB[i]);	\
	POPCOUNT2(this->helper_simd.counters[2], __intermediate);	\
	__intermediate  = PHASED_ALTREF(vectorA[i], vectorB[i]);	\
	POPCOUNT2(this->helper_simd.counters[1], __intermediate);	\
	i += 1;														\
}

#define ITER {													\
	__intermediate  = PHASED_ALTALT(vectorA[i], vectorB[i]);	\
	POPCOUNT2(this->helper_simd.counters[3], __intermediate);	\
	ITER_SHORT													\
}

	for( ; i < frontBonus; ) 					  	ITER_SHORT
	for( ; i < this->vectorCycles - tailBonus; )  	ITER
	for( ; i < this->vectorCycles - tailSmallest; ) ITER_SHORT

#undef ITER
#undef ITER_SHORT
	U32 k = this->byteAlignedEnd;
#else
	U32 k = 0;
#endif


#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k+8 < this->byte_width; k += 8){
#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
		for(U32 l = 0; l < 8; ++l){
			this->helper_simd.scalarA[l] = (arrayA[k+l] & arrayB[k+l]);
			this->helper_simd.scalarC[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayA[k+l]);
			this->helper_simd.scalarD[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayB[k+l]);
		}

		this->helper_simd.counters[2] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(this->helper_simd.scalarC));
		this->helper_simd.counters[1] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(this->helper_simd.scalarD));
		this->helper_simd.counters[3] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(this->helper_simd.scalarA));
	}

#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k < this->byte_width; ++k){
		this->helper_simd.counters[2] += POPCOUNT_ITER((arrayA[k] ^ arrayB[k]) & arrayA[k]);
		this->helper_simd.counters[1] += POPCOUNT_ITER((arrayA[k] ^ arrayB[k]) & arrayB[k]);
		this->helper_simd.counters[3] += POPCOUNT_ITER(arrayA[k] & arrayB[k]);
	}

	this->helper[1] = this->helper_simd.counters[1];
	this->helper[4] = this->helper_simd.counters[2];
	this->helper[5] = this->helper_simd.counters[3];
	this->helper[0] = this->samples*2 - (this->helper[1] + this->helper[4] + this->helper[5]);
	this->helper_simd.counters[1] = 0;
	this->helper_simd.counters[2] = 0;
	this->helper_simd.counters[3] = 0;

#if __DEBUG_MODE__ == 4
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cerr << "V\t" << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << '\t' << this->helper[0] << '\t' << this->helper[1] << '\t' << this->helper[4] << '\t' << this->helper[5] << "\t" << this->helper[0]+this->helper[1]+this->helper[4]+this->helper[5] << "\t" << ticks_per_iter.count() << std::endl;
#elif __DEBUG_MODE__ == 6
	std::cerr << "PVC\t"
			  << this->helper[0] << '\t'
			  << this->helper[1] << '\t'
			  << this->helper[4] << '\t'
			  << this->helper[5] << std::endl;
#endif

	this->setFLAGs(a, b);

	return(this->CalculateLDPhasedMath());
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDPhased(const controller_type& a, const controller_type& b){
	if(a.currentMeta().MAF == 0 || b.currentMeta().MAF == 0)
		return false;

	this->helper.resetPhased();
#if __DEBUG_MODE__ == 4
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	T currentLengthA = a[0].runs;
	T currentLengthB = b[0].runs;

	BYTE currentMixL = ((  a[0].alleleA & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)
						^ (b[0].alleleA & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1));
	BYTE currentMixR = ((  a[0].alleleB & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)
						^ (b[0].alleleB & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1));

	U32 pointerA = 0;
	U32 pointerB = 0;
	U32 add;

#if __DEBUG_MODE__ == 5
	U64 iterations = 0;
#endif

	while(true){
		if(currentLengthA > currentLengthB){ // If processed run length A > processed run length B
			currentLengthA -= currentLengthB;
			add = currentLengthB;
			++pointerB;
			currentLengthB = b[pointerB].runs;
		} else if(currentLengthA < currentLengthB) { // If processed run length A < processed run length B
			currentLengthB -= currentLengthA;
			add = currentLengthA;
			++pointerA;
			currentLengthA = a[pointerA].runs;
		} else { // If processed run length A == processed run length B
			add = currentLengthB;
			++pointerA;
			++pointerB;
			currentLengthA = a[pointerA].runs;
			currentLengthB = b[pointerB].runs;
		}
		this->helper[currentMixL] += add;
		this->helper[currentMixR] += add;

		// Exit condition
		if(pointerA == a.meta[a.metaPointer].runs || pointerB == b.meta[b.metaPointer].runs){
			if(pointerA != a.meta[a.metaPointer].runs || pointerB != b.meta[b.metaPointer].runs){
				std::cerr << Tomahawk::Helpers::timestamp("FATAL") << "Failed to exit equally!\n" << pointerA << "/" << a.meta[a.metaPointer].runs << " and " << pointerB << "/" << b.meta[b.metaPointer].runs << std::endl;
				exit(1);
			}
			break;
		}

		// Update mixing value
		currentMixL = ((  a[pointerA].alleleA & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)
					   ^ (b[pointerB].alleleA & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1));
		currentMixR = ((  a[pointerA].alleleB & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)
					   ^ (b[pointerB].alleleB & ((1 << Constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1));

	#if __DEBUG_MODE__ == 5
		++iterations;
	#endif
	}

#if __DEBUG_MODE__ == 4
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cerr << "T\t" << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << '\t' << this->helper[0] << '\t' << this->helper[1] << '\t' << this->helper[4] << '\t' << this->helper[5] << '\t' << this->helper[0]+this->helper[1]+this->helper[4]+this->helper[5] << '\t' << ticks_per_iter.count() << std::endl;
#elif __DEBUG_MODE__ == 6
	std::cerr << "P\t"
			  << this->helper[0] << '\t'
			  << this->helper[1] << '\t'
			  << this->helper[4] << '\t'
			  << this->helper[5] << std::endl;
#endif

#if __DEBUG_MODE__ == 5
	std::cerr << a.currentMeta().runs << '\t' << b.currentMeta().runs << '\t' << iterations << std::endl;
#endif

	this->setFLAGs(a, b);

	return(this->CalculateLDPhasedMath());
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDPhasedMath(void){
	// Trigger phased flag
	this->helper.setPhased();

	// Total amount of non-missing alleles
	this->helper.totalAlleleCounts = this->helper[0] + this->helper[1] + this->helper[4] + this->helper[5];

	// All values are missing
	if(this->helper.totalAlleleCounts < MINIMUM_ALLOWED_ALLELES){
		++this->insufficent_alleles;
		return false;
	}
	++this->possible;

	// Find largest
	if(this->helper.countAlternatives() < this->parameters.minimum_alleles){
		//std::cerr << "insufficient: " << max << '\t' << max2 << '\t' <<  this->parameters.minimum_alleles << std::endl;
		return false;
	}

	// Haplotype frequencies
	this->helper.haplotypeCounts[0] = (helper[0] + helper[1]) / this->helper.totalAlleleCounts;
	this->helper.haplotypeCounts[1] = (helper[4] + helper[5]) / this->helper.totalAlleleCounts;
	this->helper.haplotypeCounts[2] = (helper[0] + helper[4]) / this->helper.totalAlleleCounts;
	this->helper.haplotypeCounts[3] = (helper[1] + helper[5]) / this->helper.totalAlleleCounts;

	this->helper.D = this->helper[0]/this->helper.totalAlleleCounts * this->helper[5]/this->helper.totalAlleleCounts - this->helper[1]/this->helper.totalAlleleCounts * this->helper[4]/this->helper.totalAlleleCounts;
	this->helper.R2 = this->helper.D*this->helper.D / (((this->helper.haplotypeCounts[0] > 0 ? this->helper.haplotypeCounts[0] : 1)  * (this->helper.haplotypeCounts[1] > 0 ? this->helper.haplotypeCounts[1] : 1) * (this->helper.haplotypeCounts[2] > 0 ? this->helper.haplotypeCounts[2] : 1) * (this->helper.haplotypeCounts[3] > 0 ? this->helper.haplotypeCounts[3] : 1)));

	if(this->helper.R2 >= this->parameters.R2_min && this->helper.R2 <= this->parameters.R2_max){
		if(this->helper.D >= 0){
			this->helper.Dmax = this->helper.haplotypeCounts[0]*this->helper.haplotypeCounts[3] < this->helper.haplotypeCounts[1]*this->helper.haplotypeCounts[2]
					? this->helper.haplotypeCounts[0]*this->helper.haplotypeCounts[3]
					: this->helper.haplotypeCounts[1]*this->helper.haplotypeCounts[2];
		} else {
			this->helper.Dmax = this->helper.haplotypeCounts[0]*this->helper.haplotypeCounts[2] < this->helper.haplotypeCounts[1]*this->helper.haplotypeCounts[3]
					? -this->helper.haplotypeCounts[0]*this->helper.haplotypeCounts[2]
					: -this->helper.haplotypeCounts[1]*this->helper.haplotypeCounts[3];
		}
		this->helper.Dprime = this->helper.D / this->helper.Dmax;

		// Calculate P: Fisher's exact test
		if(this->helper.D < 0)
			this->helper.P = this->fisherController.fisherTestLess(this->helper[0],this->helper[1],this->helper[4],this->helper[5]);
		else
			this->helper.P = this->fisherController.fisherTestGreater(this->helper[0],this->helper[1],this->helper[4],this->helper[5]);

		// Toggle Fisher's exact test FLAG
		this->helper.setFisherTest();

		// Fisher's exact test P value filter
		if(this->helper.P > this->parameters.P_threshold)
			return false;

		// Calculate Chi-Sq CV from 2x2 contingency table
		this->helper.chiSqModel = 0;
		this->helper.chiSqFisher = this->fisherController.chiSquaredTest(this->helper[0],this->helper[1],this->helper[4],this->helper[5]); // Todo; fix

		return true;
	}
	return false;
}

// Execute diagonal working order
template <class T>
bool TomahawkCalculateSlave<T>::DiagonalWorkOrder(const order_type& order){
	for(U32 i = order.fromRow; i < order.toRow; ++i){
		controller_type block1(this->manager[i]);

		for(U32 j = i; j < order.toColumn; ++j){
			//std::cerr << Helpers::timestamp("DEBUG", "DIAG") << i << '/' << j << '\t' << order << std::endl;
			if(i == j)
				this->CompareBlocks(block1);
			else {
				controller_type block2(this->manager[j]);
				this->CompareBlocks(block1, block2);
			}
		}
	}
	return true;
}

// Execute square working order
template <class T>
bool TomahawkCalculateSlave<T>::SquareWorkOrder(const order_type& order){
	if(order.staggered)
		return(this->DiagonalWorkOrder(order));

	for(U32 i = order.fromRow; i < order.toRow; ++i){
		controller_type block1(this->manager[i]);

		for(U32 j = order.fromColumn; j < order.toColumn; ++j){
			//std::cerr << Helpers::timestamp("DEBUG", "SQUARE") << i << '/' << j << '\t' << order << std::endl;
			if(i == j)
				this->CompareBlocks(block1);
			else {
				controller_type block2(this->manager[j]);
				this->CompareBlocks(block1, block2);
			}
		}
	}
	return true;
}

template <class T>
bool TomahawkCalculateSlave<T>::Calculate(void){
	// If there is no data
	if(this->manager.size() == 0){
		std::cerr << Helpers::timestamp("ERROR", "CONTROLLER") << "There is no data..." << std::endl;
		return false;
	}

	// If there is no data
	if(this->orders.size() == 0){
		std::cerr << Helpers::timestamp("WARNING", "SLAVE") << "No data was given to this worker. Balancing incomplete..." << std::endl;
		return false;
	}

	// Foreach work order
	for(U32 i = 0; i < this->orders.size(); ++i){
		// If the order includes part of the diagonal
		if(!this->SquareWorkOrder(this->orders[i]))
			return false;
	}

	// Finish the output manager
	this->output_manager.Finalise();

	return true;
}

template <class T>
void TomahawkCalculateSlave<T>::CompareBlocksFunction(const controller_type& block1, const controller_type block2){
	// Method 1: Input-specified (default)
	// Method 2: Phased Vectorized No-Missing
	// Method 3:
	// Method 4: Unphased regular and unphased vectorized
	// Method 5:
	// Method 6: All algorithms comparison
	#if __DEBUG_MODE__ == 1 // 1 = No debug mode
	// Ignore when one or both is invariant
	if(block1.currentMeta().MAF == 0 || block2.currentMeta().MAF == 0 || block1.currentMeta().runs == 1 || block2.currentMeta().runs == 1){
		//std::cerr << "invariant" << std::endl;
		return;
	}

	if(block1.currentMeta().phased == 1 && block2.currentMeta().phased == 1){
		if(block1.currentMeta().MAF+block2.currentMeta().MAF <= 0.004792332){
			if(this->CalculateLDPhased(block1, block2))
				this->output_manager.Add(block1, block2, this->helper);
		} else {
			if(this->CalculateLDPhasedVectorized(block1, block2))
				this->output_manager.Add(block1, block2, this->helper);
		}
	} else {
		if(block1.currentMeta().MAF+block2.currentMeta().MAF <= 0.009784345){
			if(this->CalculateLDUnphased(block1, block2))
				this->output_manager.Add(block1, block2, this->helper);
		} else {
			if(this->CalculateLDUnphasedVectorized(block1, block2))
				this->output_manager.Add(block1, block2, this->helper);
		}
	}

	#elif __DEBUG_MODE__ == 2
	if(this->CalculateLDPhasedVectorizedNoMissing(block1, block2)){
		this->output_manager.Add(block1, block2, this->helper);
	}
	#elif __DEBUG_MODE__ == 3
		return;

	#elif __DEBUG_MODE__ == 4 | __DEBUG_MODE__ == 5
	// DEBUG
	this->CalculateLDUnphased(block1, block2);
	this->CalculateLDUnphasedVectorized(block1, block2);
	#elif __DEBUG_MODE__ == 6
	// Todo: add comprehensive
	// Every method after one another
	this->CalculateLDPhased(block1, block2);
	this->CalculateLDPhasedVectorized(block1, block2);
	this->CalculateLDPhasedVectorizedNoMissing(block1, block2);
	this->CalculateLDUnphased(block1, block2);
	this->CalculateLDUnphasedVectorized(block1, block2);
	this->CalculateLDUnphasedVectorizedNoMissing(block1, block2);
	#endif
}

template <class T>
void TomahawkCalculateSlave<T>::CompareBlocksFunctionForcedPhased(const controller_type& block1, const controller_type block2){
	// Ignore when one or both is invariant
	if(block1.currentMeta().MAF == 0 || block2.currentMeta().MAF == 0 || block1.currentMeta().runs == 1 || block2.currentMeta().runs == 1){
		//std::cerr << "invariant" << std::endl;
		return;
	}

	if(block1.currentMeta().MAF+block2.currentMeta().MAF <= 0.004792332){
		if(this->CalculateLDPhased(block1, block2))
			this->output_manager.Add(block1, block2, this->helper);
	} else {
		if(this->CalculateLDPhasedVectorized(block1, block2))
			this->output_manager.Add(block1, block2, this->helper);
	}

}

template <class T>
void TomahawkCalculateSlave<T>::CompareBlocksFunctionForcedUnphased(const controller_type& block1, const controller_type block2){
	// Ignore when one or both is invariant
	if(block1.currentMeta().MAF == 0 || block2.currentMeta().MAF == 0 || block1.currentMeta().runs == 1 || block2.currentMeta().runs == 1){
		//std::cerr << "invariant" << std::endl;
		return;
	}

	if(block1.currentMeta().MAF+block2.currentMeta().MAF <= 0.009784345){
		if(this->CalculateLDUnphased(block1, block2))
			this->output_manager.Add(block1, block2, this->helper);
	} else {
		if(this->CalculateLDUnphasedVectorized(block1, block2))
			this->output_manager.Add(block1, block2, this->helper);
	}
}

// Within-block comparisons
template <class T>
bool TomahawkCalculateSlave<T>::CompareBlocks(controller_type& block1){
	//std::cerr << Helpers::timestamp("DEBUG", "DIAG-INTERNAL") << *block1.support << '\t' << (block1.size()*block1.size()-block1.size())/2 << std::endl;
	block1.reset(); // make sure it is reset
	controller_type block2(block1);

	for(U32 i = 0; i < block1.size(); ++i){
 		block2 = block1;
		++block2; // block2 starts at relative +1
		for(U32 j = i + 1; j < block2.size(); ++j){
			(this->*phase_function_across)(block1, block2);
			++block2;
		}

		// Update progress
		this->progress(block1.size() - (i + 1), this->output_manager.GetProgressCounts());
		this->output_manager.ResetProgress();
		++block1;
	}
	return true;
}

// Across block comparisons
template <class T>
bool TomahawkCalculateSlave<T>::CompareBlocks(controller_type& block1, controller_type block2){
	// Reset
	// Make sure pointers are the beginning
	block1.reset();
	block2.reset();

	// Cycle over block 1 and block 2
	for(U32 i = 0; i < block1.size(); ++i){
		for(U32 j = 0; j < block2.size(); ++j){
			(this->*phase_function_across)(block1, block2);
			++block2;
		}

		// Update progress
		this->progress(block2.size(), this->output_manager.GetProgressCounts());
		this->output_manager.ResetProgress();

		// Reset position in block2 and increment position in block1
		block2.reset();
		++block1;
	}
	return true;
}

}

#endif /* TOMAHAWK_TOMAHAWKCALCULATESLAVE_H_ */
