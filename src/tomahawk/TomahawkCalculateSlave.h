#ifndef TOMAHAWK_TOMAHAWKCALCULATESLAVE_H_
#define TOMAHAWK_TOMAHAWKCALCULATESLAVE_H_

#include <cmath>
#include <thread>
#include <cassert>
#include <cmath>

#include "../support/simd_definitions.h"
#include "../algorithm/spinlock.h"
#include "TomahawkOutput/TomahawkOutputLD.h"
#include "../interface/ProgressBar.h"
#include "TomahawkCalcParameters.h"
#include "../math/FisherMath.h"
#include "../algorithm/LoadBalancerBlock.h"
#include "../algorithm/GenotypeBitPacker.h"
#include "TomahawkSlaveSIMDHelper.h"
#include "../io/BasicWriters.h"
#include "TomahawkIteratorManager.h"
#include "TomahawkOutput/TomahawkOutputManager.h"

// Method 1: None: Input-specified (default)
// Method 2: Phased Vectorized No-Missing
// Method 3: Count comparisons for A1
// Method 4: Unphased A1 and unphased A2
// Method 5: Phased A1 and Phased A2
// Method 6: All algorithms comparison (debug)
// Method 7: All algorithms run-time output (debug)
#define SLAVE_DEBUG_MODE	1

namespace Tomahawk{

#if SLAVE_DEBUG_MODE == 7
#pragma pack(1)
struct __costHelper{
	// RLE_A and RLE_B for A1
	float RLE_A, RLE_B;
	// cycles A1_P and A1_U
	U32 cycles_A1_P, cycles_A1_U;
	// front, bonus, mid, bonus, tail, linear cycles (for A2P_NM and A2U_NM)
	U32 skip_front, bonus_front, mid_iterations, bonus_tail, skip_tail;
	// cycles A2_1-4
	U32 cycles_A2_P, cycles_A2_U;
	U32 cycles_A2_PM, cycles_A2_UM;

	inline void reset(void){ memset(this, 0, sizeof(__costHelper)); }

	friend void operator<<(IO::BasicBuffer& buffer, const __costHelper& entry){
		buffer.Add(reinterpret_cast<const char*>(&entry), sizeof(__costHelper));
	}
};
#elif SLAVE_DEBUG_MODE == 6
/*
 This supportive structure is only used internally for
 directly comparing the output values of the two primary
 algorithms A1 and A2 and their respective variations
 */
struct __methodCompare{
	typedef __methodCompare self_type;
	typedef Tomahawk::Support::TomahawkOutputLD helper_type;

	__methodCompare(){}
	~__methodCompare(){}

	float phased[3][4];  // Phased A1, phased no missing A2, phased missing A2
	float unphased[3][9]; // Unphased A1, unphased no missing A2, unphased missing A2

	friend std::ostream& operator<<(std::ostream& os, const self_type& m){
		// P, PV, PVM, U, UV, UVM
		os << "P\t" << m.phased[0][0] << '\t' << m.phased[0][1] << '\t' << m.phased[0][2] << '\t' << m.phased[0][3] << std::endl;
		os << "PV\t" << m.phased[1][0] << '\t' << m.phased[1][1] << '\t' << m.phased[1][2] << '\t' << m.phased[1][3] << std::endl;
		os << "PVM\t" << m.phased[2][0] << '\t' << m.phased[2][1] << '\t' << m.phased[2][2] << '\t' << m.phased[2][3] << std::endl;
		os << "U\t";
		for(U32 i = 0; i < 8; ++i) os << m.unphased[0][i] << '\t'; os << m.unphased[0][8] << std::endl;;
		os << "UV\t";
		for(U32 i = 0; i < 8; ++i) os << m.unphased[1][i] << '\t'; os << m.unphased[1][8] << std::endl;;
		os << "UVM\t";
		for(U32 i = 0; i < 8; ++i) os << m.unphased[2][i] << '\t'; os << m.unphased[2][8] << std::endl;;

		return(os);
	}

	void addPhased(const U32 p, const helper_type& helper){
		this->phased[p][0] = helper[0];
		this->phased[p][1] = helper[1];
		this->phased[p][2] = helper[4];
		this->phased[p][3] = helper[5];
	}

	void addUnphased(const U32 p, const helper_type& helper){
		this->unphased[p][0] = helper[0];
		this->unphased[p][1] = helper[1] + helper[4];
		this->unphased[p][2] = helper[5];
		this->unphased[p][3] = helper[16] + helper[64];
		this->unphased[p][4] = helper[17] + helper[20] + helper[65] + helper[68];
		this->unphased[p][5] = helper[21] + helper[69];
		this->unphased[p][6] = helper[80];
		this->unphased[p][7] = helper[81] + helper[84];
		this->unphased[p][8] = helper[85];
	}

	// Check to make sure all algorithms and variations
	// produce the correct output values
	bool validate(void) const{
		for(U32 i = 0; i < 4; ++i){
			if(this->phased[0][i] != this->phased[1][i] ||
			   this->phased[0][i] != this->phased[2][i] ||
			   this->phased[1][i] != this->phased[2][i])
			{
				std::cerr << Helpers::timestamp("ERROR", "VALIDATION") << "Phased failure: " << this->phased[0][i]
						  << '\t' << this->phased[1][i] << '\t' << this->phased[2][i] << std::endl;
				return false;
			}
		}

		for(U32 i = 0; i < 9; ++i){
			if(this->unphased[0][i] != this->unphased[1][i] ||
			   this->unphased[0][i] != this->unphased[2][i] ||
			   this->unphased[1][i] != this->unphased[2][i])
			{
				std::cerr << Helpers::timestamp("ERROR", "VALIDATION") << "Phased failure: " << this->unphased[0][i]
						  << '\t' << this->unphased[1][i] << '\t' << this->unphased[2][i] << std::endl;
				return false;
			}
		}

		return true;
	}
};

#endif

// Parameter thresholds for FLAGs
#define LOW_MAF_THRESHOLD		0.01
#define LOW_HWE_THRESHOLD		1e-6
#define LONG_RANGE_THRESHOLD	500e3
#define MINIMUM_ALLOWED_ALLELES	5		// Minimum number of alleles required for math to work out in the unphased case
#define CHI_SQ_MAX_CV			1300

// SIMD trigger
#if SIMD_AVAILABLE == 1

#ifdef _popcnt64
#define POPCOUNT_ITER	_popcnt64
#else
#define POPCOUNT_ITER	__builtin_popcountll
#endif

#define UNPHASED_UPPER_MASK	170  // 10101010b
#define UNPHASED_LOWER_MASK	85   // 01010101b
#define FILTER_UNPHASED_BYTE(A, B) (((((A & UNPHASED_UPPER_MASK) | (B & UNPHASED_LOWER_MASK)) & UNPHASED_LOWER_MASK) << 1) & A)
#define FILTER_UNPHASED_BYTE_PAIR(A, B, C, D) ((FILTER_UNPHASED_BYTE(A, B) >> 1) | FILTER_UNPHASED_BYTE(C, D))
#define FILTER_UNPHASED_BYTE_SPECIAL(A) (((A >> 1) & A) & UNPHASED_LOWER_MASK)

#if SIMD_VERSION == 6 // AVX-512: UNTESTED
const VECTOR_TYPE ONE_MASK = _mm512_set1_epi8(255); // 11111111b
const VECTOR_TYPE maskUnphasedHigh = _mm512_set1_epi8(UNPHASED_UPPER_MASK);	// 10101010b
const VECTOR_TYPE maskUnphasedLow  = _mm512_set1_epi8(UNPHASED_LOWER_MASK);	// 01010101b

#define PHASED_ALTALT(A,B) _mm512_and_si512(A, B)
#define PHASED_REFREF(A,B) _mm512_and_si512(_mm512_xor_si512(A, ONE_MASK), _mm512_xor_si512(B, ONE_MASK))
#define PHASED_ALTREF(A,B) _mm512_and_si512(_mm512_xor_si512(A, B), B)
#define PHASED_REFALT(A,B) _mm512_and_si512(_mm512_xor_si512(A, B), A)
#define PHASED_ALTALT_MASK(A,B,M) _mm512_and_si512(PHASED_ALTALT(A, B), M)
#define PHASED_REFREF_MASK(A,B,M) _mm512_and_si512(PHASED_REFREF(A, B), M)
#define PHASED_ALTREF_MASK(A,B,M) _mm512_and_si512(PHASED_ALTREF(A, B), M)
#define PHASED_REFALT_MASK(A,B,M) _mm512_and_si512(PHASED_REFALT(A, B), M)
#define MASK_MERGE(A,B) _mm512_xor_si512(_mm512_or_si512(A, B), ONE_MASK)

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

#if SIMD_VERSION >= 3
#define POPCOUNT(A, B) {								\
	A += POPCOUNT_ITER(_mm_extract_epi64(B, 0));		\
	A += POPCOUNT_ITER(_mm_extract_epi64(B, 1));		\
}
#else
#define POPCOUNT(A, B) { \
	U64 temp = _mm_extract_epi16(B, 0) << 6 | _mm_extract_epi16(B, 1) << 4 | _mm_extract_epi16(B, 2) << 2 | _mm_extract_epi16(B, 3); \
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

template <class T>
class TomahawkCalculateSlave{
	//Basic typedefs
	typedef TomahawkCalculateSlave<T> self_type;
	typedef TomahawkBlockManager<const T> manager_type;
	typedef TomahawkIterator<const T> controller_type;
	typedef const Support::TomahawkEntryMeta<const T> meta_type;
	typedef const Support::TomahawkRun<const T> run_type;
	typedef Totempole::TotempoleEntry totempole_entry_type;
	typedef IO::TomahawkOutputManager<T> output_manager_type;
	typedef Support::TomahawkOutputLD helper_type;
	typedef TomahawkBlockPackedPair<> simd_pair;

	// Work orders
	typedef Tomahawk::LoadBalancerBlock order_type;
	typedef std::vector<order_type> work_order;

	// Function pointers
	typedef void (self_type::*phaseFunction)(const controller_type& block1, const controller_type block2);

public:
	TomahawkCalculateSlave(const manager_type& manager,
		output_manager_type& writer,
		Interface::ProgressBar& progress,
		const TomahawkCalcParameters& parameters,
		const work_order& orders);

	~TomahawkCalculateSlave();

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
	output_manager_type& getOutputManager(void){ return(this->output_manager); }

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
	Algorithm::FisherMath fisherController;
	const manager_type& manager;

	// thread
	std::thread thread;

	// writer manager
	output_manager_type output_manager; // each thread has their own output manager

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
TomahawkCalculateSlave<T>::TomahawkCalculateSlave(const manager_type& manager,
		output_manager_type& writer,
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
	output_manager(writer),
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
	this->output_manager += other.output_manager;
	return(*this);
}

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
	if(mA.MGF < LOW_MAF_THRESHOLD)
		this->helper.setLowMAFA();

	if(mB.MGF < LOW_MAF_THRESHOLD)
		this->helper.setLowMAFB();

	if(mA.HWE_P < LOW_HWE_THRESHOLD)
		this->helper.setFailedHWEA();

	if(mB.HWE_P < LOW_HWE_THRESHOLD)
		this->helper.setFailedHWEB();

	if(mA.missing || mB.missing)
		this->helper.setHasMissingValues();
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDUnphased(const controller_type& a, const controller_type& b){
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

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	while(true){
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

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\t';
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

	// For writing output
	this->helper[0] = this->helper.haplotypeCounts[0] * 2*this->helper.totalAlleleCounts;
	this->helper[1] = this->helper.haplotypeCounts[1] * 2*this->helper.totalAlleleCounts;
	this->helper[4] = this->helper.haplotypeCounts[2] * 2*this->helper.totalAlleleCounts;
	this->helper[5] = this->helper.haplotypeCounts[3] * 2*this->helper.totalAlleleCounts;

	this->helper.D = this->helper.haplotypeCounts[0] * this->helper.haplotypeCounts[3] - this->helper.haplotypeCounts[1] * this->helper.haplotypeCounts[2];
	this->helper.R2 = this->helper.D*this->helper.D / (p * (1 - p) * q * (1 - q));

	if(this->helper.countAlternatives() < this->parameters.minimum_alleles)
		return false;

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

		 if(this->helper.D < 0)
			this->helper.P = this->fisherController.fisherTestLess(round(this->helper[0]),round(this->helper[1]),round(this->helper[4]),round(this->helper[5]));
		else
			this->helper.P = this->fisherController.fisherTestGreater(round(this->helper[0]),round(this->helper[1]),round(this->helper[4]),round(this->helper[5]));

		if(this->helper[0] < 1 || this->helper[1] < 1 || this->helper[4] < 1 || this->helper[5] < 1)
			this->helper.setIncomplete();

		// Fisher's exact test P value filter
		if(this->helper.P > this->parameters.P_threshold)
			return false;

		this->helper.chiSqFisher = this->fisherController.chiSquaredTest(this->helper[0],this->helper[1],this->helper[4],this->helper[5]);

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
		const float p0 = 2*this->helper[0] + this->helper[1]  + this->helper[4]    + this->helper[16] + this->helper[64];
		const float q0 = this->helper[16]  + this->helper[64] + 2*this->helper[80] + this->helper[81] + this->helper[84];
		const float p1 = this->helper[1]   + this->helper[4]  + 2*this->helper[5]  + this->helper[21] + this->helper[69];
		const float q1 = this->helper[21]  + this->helper[69] + this->helper[81]   + this->helper[84] + 2*this->helper[85];

		this->helper[0] = p0;
		this->helper[1] = p1;
		this->helper[4] = q0;
		this->helper[5] = q1;

		// Update counter
		++this->no_uncertainty;

		// Reset
		this->helper.chiSqModel = 0;

		// Use standard math
		return(this->CalculateLDPhasedMath());
	}

	const double p = ((this->helper[0] + this->helper[1]  + this->helper[4]  + this->helper[5])*2.0
				   + (this->helper[16] + this->helper[17] + this->helper[20] + this->helper[21] + this->helper[64] + this->helper[65] + this->helper[68] + this->helper[69]))
				   / (2.0 * this->helper.totalAlleleCounts);
	const double q = ((this->helper[0] + this->helper[16] + this->helper[64] + this->helper[80])*2.0
				   + (this->helper[1]  + this->helper[4]  + this->helper[17] + this->helper[20] + this->helper[65] + this->helper[68] + this->helper[81] + this->helper[84]))
				   / (2.0 * this->helper.totalAlleleCounts);
	const double n11 = 2.0* this->helper[0] + this->helper[1] + this->helper[4] + this->helper[16] + this->helper[64];

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
#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
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
	POPCOUNT(this->helper_simd.counters[0], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);	\
	POPCOUNT(this->helper_simd.counters[1], __intermediate);				\
	__intermediate = FILTER_UNPHASED(altref, altref);						\
	POPCOUNT(this->helper_simd.counters[2], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);	\
	POPCOUNT(this->helper_simd.counters[3], __intermediate);				\
	__intermediate = FILTER_UNPHASED(refalt,refalt);						\
	POPCOUNT(this->helper_simd.counters[6], __intermediate);				\
}

#define ITER_LONG {															\
	ITER_SHORT																\
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref);	\
	POPCOUNT(this->helper_simd.counters[5], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt);	\
	POPCOUNT(this->helper_simd.counters[7], __intermediate);				\
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt);						\
	POPCOUNT(this->helper_simd.counters[8], __intermediate);				\
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

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << ticks_per_iter.count() << '\n';
#endif

	this->setFLAGs(a, b);
	return(this->CalculateLDUnphasedMath());
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDUnphasedVectorized(const controller_type& a, const controller_type& b){
	#if SLAVE_DEBUG_MODE < 6
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
#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
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
	POPCOUNT(this->helper_simd.counters[0], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);	\
	POPCOUNT(this->helper_simd.counters[1], __intermediate);				\
	__intermediate = FILTER_UNPHASED(altref, altref); 						\
	POPCOUNT(this->helper_simd.counters[2], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);	\
	POPCOUNT(this->helper_simd.counters[3], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref, altalt, altalt, refref);	\
	POPCOUNT(this->helper_simd.counters[4], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altref, altref, refalt);	\
	POPCOUNT(this->helper_simd.counters[4], __intermediate);				\
	__intermediate = FILTER_UNPHASED(refalt,refalt);						\
	POPCOUNT(this->helper_simd.counters[6], __intermediate);				\
}

#define ITER_LONG {															\
	ITER_SHORT																\
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref); 	\
	POPCOUNT(this->helper_simd.counters[5], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt); 	\
	POPCOUNT(this->helper_simd.counters[7], __intermediate);				\
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt); 						\
	POPCOUNT(this->helper_simd.counters[8], __intermediate);				\
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

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << ticks_per_iter.count() << '\n';
#endif

	this->setFLAGs(a, b);
	return(this->CalculateLDUnphasedMath());
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDPhasedVectorized(const controller_type& a, const controller_type& b){
	#if SLAVE_DEBUG_MODE < 6
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
#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_SHORT {														\
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);					\
	__intermediate  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(this->helper_simd.counters[0], __intermediate);				\
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(this->helper_simd.counters[1], __intermediate);				\
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(this->helper_simd.counters[2], __intermediate);				\
	i += 1;																	\
}

#define ITER {																\
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);					\
	__intermediate  = PHASED_ALTALT_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(this->helper_simd.counters[3], __intermediate);				\
	__intermediate  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(this->helper_simd.counters[0], __intermediate);				\
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(this->helper_simd.counters[1], __intermediate);				\
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(this->helper_simd.counters[2], __intermediate);				\
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
	this->helper[0] = (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 + this->helper_simd.counters[0] - this->phased_unbalanced_adjustment;


#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << "V\t" << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif

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
	const U32 tailSmallest  = datA.tailZero  < datB.tailZero  ? datA.tailZero  : datB.tailZero;
	U32 i = frontSmallest;
	const U32 frontBonus = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus  = datA.tailZero  != tailSmallest  ? datA.tailZero  : datB.tailZero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	//VECTOR_TYPE altalt, altref, refalt;
	VECTOR_TYPE __intermediate;

// Debug timings
#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_SHORT {											\
	__intermediate  = PHASED_REFALT(vectorA[i], vectorB[i]);	\
	POPCOUNT(this->helper_simd.counters[2], __intermediate);	\
	__intermediate  = PHASED_ALTREF(vectorA[i], vectorB[i]);	\
	POPCOUNT(this->helper_simd.counters[1], __intermediate);	\
	i += 1;														\
}

#define ITER {													\
	__intermediate  = PHASED_ALTALT(vectorA[i], vectorB[i]);	\
	POPCOUNT(this->helper_simd.counters[3], __intermediate);	\
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

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << "V\t" << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << "\t" << ticks_per_iter.count() << '\n';
#endif

	this->setFLAGs(a, b);

	return(this->CalculateLDPhasedMath());
}

template <class T>
bool TomahawkCalculateSlave<T>::CalculateLDPhased(const controller_type& a, const controller_type& b){
	this->helper.resetPhased();
#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
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

#if SLAVE_DEBUG_MODE == 3
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

	#if SLAVE_DEBUG_MODE == 3
		++iterations;
	#endif
	}

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << a.currentMeta().MAF*this->samples + b.currentMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif

#if SLAVE_DEBUG_MODE == 3
	std::cout << a.currentMeta().runs << '\t' << b.currentMeta().runs << '\t' << iterations << std::endl;
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

	// Filter by minor haplotype frequency
	if(this->helper.countAlternatives() < this->parameters.minimum_alleles)
		return false;


	// Haplotype frequencies
	this->helper.haplotypeCounts[0] = (this->helper[0] + this->helper[1]) / this->helper.totalAlleleCounts;
	this->helper.haplotypeCounts[1] = (this->helper[4] + this->helper[5]) / this->helper.totalAlleleCounts;
	this->helper.haplotypeCounts[2] = (this->helper[0] + this->helper[4]) / this->helper.totalAlleleCounts;
	this->helper.haplotypeCounts[3] = (this->helper[1] + this->helper[5]) / this->helper.totalAlleleCounts;

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

		if(this->helper[0] == 0 || this->helper[1] == 0 || this->helper[4] == 0 || this->helper[5] == 0)
			this->helper.setIncomplete();

		// Fisher's exact test P value filter
		if(this->helper.P > this->parameters.P_threshold){
			//std::cerr << this->helper.P << '\t' << this->helper.D << '\t' << this->helper[0] << '\t' << this->helper[1] << '\t' << this->helper[4] << '\t' << this->helper[5] << std::endl;
			//if(this->helper.P > 1)
			//	exit(1);

			return false;
		}

		// Calculate Chi-Sq CV from 2x2 contingency table
		this->helper.chiSqModel = 0;
		this->helper.chiSqFisher = this->fisherController.chiSquaredTest(this->helper[0],this->helper[1],this->helper[4],this->helper[5]);

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
	this->output_manager.flushBlock();

	return true;
}

template <class T>
void TomahawkCalculateSlave<T>::CompareBlocksFunction(const controller_type& block1, const controller_type block2){
#if SLAVE_DEBUG_MODE == 1 // 1 = No debug mode
	// Ignore when one or both is invariant
	if(block1.currentMeta().runs == 1 || block2.currentMeta().runs == 1){
		//std::cerr << "invariant" << std::endl;
		return;
	}

	if(block1.currentMeta().phased == 1 && block2.currentMeta().phased == 1){
		if(block1.currentMeta().MGF + block2.currentMeta().MGF <= 0.004792332){
			if(this->CalculateLDPhased(block1, block2))
				this->output_manager.Add(block1, block2, this->helper);
		} else {
			if(this->CalculateLDPhasedVectorized(block1, block2))
				this->output_manager.Add(block1, block2, this->helper);
		}
	} else {
		if(block1.currentMeta().MGF+block2.currentMeta().MGF <= 0.009784345){
			if(this->CalculateLDUnphased(block1, block2))
				this->output_manager.Add(block1, block2, this->helper);
		} else {
			if(this->CalculateLDUnphasedVectorized(block1, block2))
				this->output_manager.Add(block1, block2, this->helper);
		}
	}

#elif SLAVE_DEBUG_MODE == 2
	if(this->CalculateLDPhasedVectorizedNoMissing(block1, block2)){
		this->output_manager.Add(block1, block2, this->helper);
	}
#elif SLAVE_DEBUG_MODE == 3
		this->CalculateLDPhased(block1, block2);

#elif SLAVE_DEBUG_MODE == 4
	// DEBUG
	this->CalculateLDUnphased(block1, block2);
	this->CalculateLDUnphasedVectorized(block1, block2);
#elif SLAVE_DEBUG_MODE == 5
	this->CalculateLDPhased(block1, block2);
	this->CalculateLDPhasedVectorized(block1, block2);
#elif SLAVE_DEBUG_MODE == 6
	// Every method after one another
	// P, PV, PVM, U, UV, UVM
	__methodCompare m;
	this->CalculateLDPhased(block1, block2);
	m.addPhased(0, this->helper);
	this->CalculateLDPhasedVectorized(block1, block2);
	m.addPhased(1, this->helper);
	this->CalculateLDPhasedVectorizedNoMissing(block1, block2);
	m.addPhased(2, this->helper);
	this->CalculateLDUnphased(block1, block2);
	m.addUnphased(0, this->helper);
	this->CalculateLDUnphasedVectorized(block1, block2);
	m.addUnphased(1, this->helper);
	this->CalculateLDUnphasedVectorizedNoMissing(block1, block2);
	m.addUnphased(2, this->helper);

	if(!m.validate()){
		std::cerr << Helpers::timestamp("ERROR", "VALIDATION") << "Failed validation..." << std::endl;
		std::cerr << m << std::endl;
		//exit(1);
	}
#endif
}

template <class T>
void TomahawkCalculateSlave<T>::CompareBlocksFunctionForcedPhased(const controller_type& block1, const controller_type block2){
	// Ignore when one or both is invariant
	if(block1.currentMeta().runs == 1 || block2.currentMeta().runs == 1){
		//std::cerr << "invariant" << std::endl;
		return;
	}

	if(block1.currentMeta().MGF + block2.currentMeta().MGF <= 0.004792332){
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
	if(block1.currentMeta().runs == 1 || block2.currentMeta().runs == 1){
		//std::cerr << "invariant" << std::endl;
		return;
	}

	if(block1.currentMeta().MGF + block2.currentMeta().MGF <= 0.009784345){
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
