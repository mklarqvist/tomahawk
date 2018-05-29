#ifndef TOMAHAWK_LD_CALCULATION_SLAVE_H_
#define TOMAHAWK_LD_CALCULATION_SLAVE_H_

#include <cmath>
#include <thread>
#include <cassert>
#include <cmath>

#include "io/basic_writers.h"
#include "support/simd_definitions.h"
#include "algorithm/spinlock.h"
#include "interface/progressbar.h"
#include "TomahawkCalcParameters.h"
#include "algorithm/load_balancer_block.h"
#include "algorithm/genotype_bitpacker.h"
#include "math/fisher_math.h"
#include "genotype_meta_container_reference.h"
#include "ld_calculation_simd_helper.h"
#include "two/output_entry_support.h"
#include "io/output_writer.h"
#include "haplotype_bitvector.h"

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
#define FILTER_UNPHASED_BYTE(A, B)            (((((A & UNPHASED_UPPER_MASK) | (B & UNPHASED_LOWER_MASK)) & UNPHASED_LOWER_MASK) << 1) & A)
#define FILTER_UNPHASED_BYTE_PAIR(A, B, C, D) ((FILTER_UNPHASED_BYTE(A, B) >> 1) | FILTER_UNPHASED_BYTE(C, D))
#define FILTER_UNPHASED_BYTE_SPECIAL(A)       (((A >> 1) & A) & UNPHASED_LOWER_MASK)

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
class LDSlave{
	typedef LDSlave<T>                                 self_type;
	typedef GenotypeMetaContainerReference<T>          manager_type;
	typedef base::GenotypeContainerReference<T>        block_type;
	typedef const support::GenotypeDiploidRun<T>       run_type;
	typedef base::GenotypeContainerRunlengthObjects<T> rle_type;

	typedef const MetaEntry                      meta_type;
	typedef totempole::IndexEntry                totempole_entry_type;
	typedef io::OutputWriter                     output_writer_type;
	typedef support::OutputEntrySupport          helper_type;
	typedef base::GenotypeBitvector<>            simd_pair;

	typedef interface::ProgressBar               progress_bar_type;
	typedef support::LDCalculationSIMDHelper<>   simd_helper_type;
	typedef TomahawkCalcParameters               parameter_type;

	// Work orders
	typedef LoadBalancerThread                   order_type;
	typedef std::vector<order_type>              work_order;

	// Function pointers
	typedef bool (self_type::*phaseFunction)(const block_type& block1, const block_type& block2);

public:
	LDSlave(const manager_type& manager,
			output_writer_type& writer,
        interface::ProgressBar& progress,
  const TomahawkCalcParameters& parameters,
              const work_order& orders);

	~LDSlave();

	LDSlave& operator=(const LDSlave& other);
	LDSlave& operator=(LDSlave&& other) noexcept;
	LDSlave& operator+=(const LDSlave& other);

	/**<
	 * Start calculations for the internal thread
	 * @return Returns a pointer to the internal thread
	 */
	inline std::thread* Start(void){
		this->thread = std::thread(&self_type::Calculate, this);
		return(&this->thread);
	}

	//inline const U64& getImpossible(void) const{ return this->impossible; }
	//inline const U64& getPossible(void) const{ return this->possible; }
	//inline const U64& getNoHets(void) const{ return this->no_uncertainty; }
	//inline const U64& getInsufficientData(void) const{ return this->insufficent_alleles; }
	//inline U64 getComparisons(void) const{ return(this->impossible + this->possible + this->insufficent_alleles); }
	inline output_writer_type& getWriter(void){ return(this->output_writer); }
	inline const output_writer_type& getWriter(void) const{ return(this->output_writer); }

private:
	/**<
	 * Main entry-point for calculating pairwise LD. This is invoked by the `Start()` function
	 * @return
	 */
	bool Calculate();

	/**<
	 * Calculate pairwise LD when the data are in the same `twk` block
	 * @param order
	 * @return
	 */
	bool DiagonalWorkOrder(const order_type& order);

	/**<
	 * Calculate pairwise LD when the data are in different `twk` blocks
	 * @param order
	 * @return
	 */
	bool SquareWorkOrder(const order_type& order);

	// Comparator functions
	bool CompareBlocks(block_type& block1);
	bool CompareBlocks(block_type& block1, block_type& block2);
	bool CompareBlocksFunction(const block_type& block1, const block_type& block2);
	bool CompareBlocksFunctionForcedPhased(const block_type& block1, const block_type& block2);
	bool CompareBlocksFunctionForcedPhasedSimple(const block_type& block1, const block_type& block2);
	bool CompareBlocksFunctionForcedPhasedSimpleSample(const block_type& block1, const block_type& block2);
	bool CompareBlocksFunctionForcedUnphased(const block_type& block1, const block_type& block2);

	// Phased functions
	bool CalculateLDPhased(const block_type& block1, const block_type& block2);
	bool CalculateLDPhasedSimple(const block_type& block1, const block_type& block2);
	bool CalculateLDPhasedSimpleSample(const block_type& block1, const block_type& block2);
	bool CalculateLDPhasedVectorized(const block_type& block1, const block_type& block2);
	bool CalculateLDPhasedVectorizedNoMissing(const block_type& block1, const block_type& block2);
	bool CalculateLDPhasedVectorizedNoMissingNoTable(const block_type& block1, const block_type& block2);

	// Unphased functions
	bool CalculateLDUnphased(const block_type& block1, const block_type& block2);
	bool CalculateLDUnphasedVectorized(const block_type& block1, const block_type& block2);
	bool CalculateLDUnphasedVectorizedNoMissing(const block_type& block1, const block_type& block2);

	// Math
	bool CalculateLDPhasedMath();
	bool CalculateLDPhasedMathSimple(const block_type& block1, const block_type& block2);
	bool CalculateLDUnphasedMath();

	// General functions
	inline const double ChiSquaredUnphasedTable(const double& target, const double& p, const double& q);
	bool ChooseF11Calculate(const double& target, const double& p, const double& q);
	inline void setFLAGs(const block_type& block1, const block_type& block2);

	// Tests

private:
	const parameter_type& parameters;

	// Counters
	U32 block_comparisons;
	U64 variant_comparisons;

	// Counters
	const U64 samples;
	//U64 impossible;
	//U64 possible;
	//U64 no_uncertainty;
	//U64 insufficent_alleles;
	//U64 false_positive;
	//U64 false_negative;

	//helper_type helper;
	simd_helper_type helper_simd;
	//Algorithm::FisherMath fisherController;
	const manager_type& manager;

	// thread
	std::thread thread;

	// writer manager
	output_writer_type output_writer; // each thread has their own output manager with its own buffer

	// progress
	progress_bar_type& progress;

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

	helper_type helper;
};

template <class T>
LDSlave<T>::LDSlave(const manager_type& manager,
		output_writer_type& writer,
		progress_bar_type& progress,
		const parameter_type& parameters,
		const work_order& orders) :
	parameters(parameters),
	block_comparisons(0),
	variant_comparisons(0),
	samples(manager.numberSamples()),
	//impossible(0),
	//possible(0),
	//no_uncertainty(0),
	//insufficent_alleles(0),
	//false_positive(0),
	//false_negative(0),
	//fisherController(1024),
	manager(manager),
	output_writer(writer),
	progress(progress),
	phase_function_across(nullptr),
	orders(orders),
	byte_width(ceil((double)samples/4)),
	byteAlignedEnd(this->byte_width/(GENOTYPE_TRIP_COUNT/4)*(GENOTYPE_TRIP_COUNT/4)),
	vectorCycles(this->byteAlignedEnd*4/GENOTYPE_TRIP_COUNT),
	phased_unbalanced_adjustment((this->samples*2)%8),
	unphased_unbalanced_adjustment(this->samples%4)
{
	// Decide block comparator function
	if(this->parameters.force == parameter_type::force_method::none)
		this->phase_function_across = &self_type::CompareBlocksFunction;
	else if(this->parameters.force == parameter_type::force_method::phasedFunction){
		if(this->parameters.fast_mode) this->phase_function_across = &self_type::CompareBlocksFunctionForcedPhasedSimple;
		else this->phase_function_across = &self_type::CompareBlocksFunctionForcedPhased;
	} else
		this->phase_function_across = &self_type::CompareBlocksFunctionForcedUnphased;
}

template <class T>
LDSlave<T>::~LDSlave(){ }

// Reduce function
template <class T>
LDSlave<T>& LDSlave<T>::operator+=(const LDSlave<T>& other){
	this->block_comparisons   += other.block_comparisons;
	this->variant_comparisons += other.variant_comparisons;
	//this->impossible          += other.impossible;
	//this->possible            += other.possible;
	//this->no_uncertainty      += other.no_uncertainty;
	//this->insufficent_alleles += other.insufficent_alleles;
	//this->false_positive    += other.false_positive;
	//this->false_negative    += other.false_negative;
	this->output_writer       += other.output_writer;
	return(*this);
}

template <class T>
inline void LDSlave<T>::setFLAGs(const block_type& block1, const block_type& block2){
	// If long range
	const meta_type& mA = block1.currentMeta();
	const meta_type& mB = block2.currentMeta();

	helper.setSameContig(block1.getTotempole().contigID == block2.getTotempole().contigID);
	helper.setLongRange(mB.position - mA.position > LONG_RANGE_THRESHOLD);
	helper.setLowAFA(mA.AF < LOW_MAF_THRESHOLD);
	helper.setLowAFB(mB.AF < LOW_MAF_THRESHOLD);
	helper.setFailedHWEA(mA.HWE_P < LOW_HWE_THRESHOLD);
	helper.setFailedHWEB(mB.HWE_P < LOW_HWE_THRESHOLD);
	helper.setHasMissingValuesA(mA.has_missing);
	helper.setHasMissingValuesB(mB.has_missing);
}

template <class T>
bool LDSlave<T>::CalculateLDUnphased(const block_type& block1, const block_type& block2){
	if(block1.currentMeta().AF == 0 || block2.currentMeta().AF == 0)
		return false;

	helper.resetUnphased();

	/*////////////
	// Calculate
	////////////*/
	const run_type* const a = block1.current();
	const run_type* const b = block2.current();
	T currentLengthA = a[0].runs;
	T currentLengthB = b[0].runs;

	BYTE currentMix = ((a[0].alleleA & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH*3)
					^ ((a[0].alleleB & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH*2)
					^ ((b[0].alleleA & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH*1)
					^ ((b[0].alleleB & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)));

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
			helper[currentMix] += currentLengthB;
			++pointerB;
			currentLengthB = b[pointerB].runs;
		}
		// If processed run length A == processed run length B
		else if(currentLengthA == currentLengthB){
			helper[currentMix] += currentLengthB;
			++pointerA;
			++pointerB;
			currentLengthA = a[pointerA].runs;
			currentLengthB = b[pointerB].runs;
		}
		// If processed run length A < processed run length B
		else {
			currentLengthB -= currentLengthA;
			helper[currentMix] += currentLengthA;
			++pointerA;
			currentLengthA = a[pointerA].runs;
		}

		// Exit condition
		if(pointerA == block1.currentMeta().runs || pointerB == block2.currentMeta().runs){
			//std::cerr << pointerA << '/' << a.meta[a.metaPointer].runs << '\t' << pointerB << '/' << b.meta[b.metaPointer].runs << std::endl;
			if(pointerA != block1.currentMeta().runs || pointerB != block2.currentMeta().runs){
				std::cerr << helpers::timestamp("FATAL") << "Failed to exit equally!\n" << pointerA << "/" << block1.currentMeta().runs << " and " << pointerB << "/" << block2.currentMeta().runs << std::endl;
				exit(1);
			}
			break;
		}

		currentMix  = ((a[pointerA].alleleA & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH*3)
					^ ((a[pointerA].alleleB & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH*2)
					^ ((b[pointerB].alleleA & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH*1)
					^ ((b[pointerB].alleleB & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)));
	}

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\t';
#endif

	this->setFLAGs(block1, block2);

	return(this->CalculateLDUnphasedMath());
}

template <class T>
inline const double LDSlave<T>::ChiSquaredUnphasedTable(const double& target, const double& p, const double& q){
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

	const double chisq1111 = e1111 > 0 ? pow(helper[0]  - e1111, 2) / e1111 : 0,
				 chisq1112 = e1112 > 0 ? pow(helper[1]  + helper[4] - e1112, 2) / e1112 : 0,
				 chisq1122 = e1122 > 0 ? pow(helper[5]  - e1122, 2) / e1122 : 0,
				 chisq1211 = e1211 > 0 ? pow(helper[16] + helper[64] - e1211, 2) / e1211 : 0,
				 chisq1212 = e1212 > 0 ? pow(helper[17] + helper[20] + helper[65] + helper[68] - e1212, 2) / e1212 : 0,
				 chisq1222 = e1222 > 0 ? pow(helper[21] + helper[69] - e1222, 2) / e1222 : 0,
				 chisq2211 = e2211 > 0 ? pow(helper[80] - e2211, 2) / e2211 : 0,
				 chisq2212 = e2212 > 0 ? pow(helper[81] + helper[84] - e2212, 2) / e2212 : 0,
				 chisq2222 = e2222 > 0 ? pow(helper[85] - e2222, 2) / e2222 : 0;

	return(chisq1111 + chisq1112 + chisq1122 + chisq1211 + chisq1212 + chisq1222 + chisq2211 + chisq2212 + chisq2222);
}

template <class T>
bool LDSlave<T>::ChooseF11Calculate(const double& target, const double& p, const double& q){
	helper.haplotypeCounts[0] = target;
	helper.haplotypeCounts[1] = p - helper.haplotypeCounts[0];
	helper.haplotypeCounts[2] = q - helper.haplotypeCounts[0];
	helper.haplotypeCounts[3] = 1 - (helper.haplotypeCounts[0] + helper.haplotypeCounts[1] + helper.haplotypeCounts[2]);

	// For writing output
	helper[0] = helper.haplotypeCounts[0] * 2*helper.totalHaplotypeCounts;
	helper[1] = helper.haplotypeCounts[1] * 2*helper.totalHaplotypeCounts;
	helper[4] = helper.haplotypeCounts[2] * 2*helper.totalHaplotypeCounts;
	helper[5] = helper.haplotypeCounts[3] * 2*helper.totalHaplotypeCounts;

	helper.D  = helper.haplotypeCounts[0] * helper.haplotypeCounts[3] - helper.haplotypeCounts[1] * helper.haplotypeCounts[2];
	helper.R2 = helper.D*helper.D / (p * (1 - p) * q * (1 - q));
	helper.R  = sqrt(helper.R2);

	if(helper.getTotalAltHaplotypeCount() < this->parameters.minimum_sum_alternative_haplotype_count)
		return false;

	if(helper.R2 < this->parameters.R2_min || helper.R2 > this->parameters.R2_max)
		return false;

	if(helper.D >= 0){
		helper.Dmax = p * (1.0-q) < q * (1.0-p)
				? p * (1.0-q)
				: q * (1.0-p);
	} else {
		helper.Dmax = p*q < (1-p)*(1-q)
				? p*q
				: (1-p)*(1-q);
	}
	helper.Dprime = helper.D / helper.Dmax;

	double left,right,both;
	kt_fisher_exact(round(helper[0]),round(helper[1]),
					round(helper[4]),round(helper[5]),
					&left,&right,&both);
	helper.P = both;

	// Fisher's exact test P value filter
	if(helper.P > this->parameters.P_threshold)
		return false;


	helper.setCompleteLD(helper[0] < 1 || helper[1] < 1 || helper[4] < 1 || helper[5] < 1);
	helper.setPerfectLD(helper.R2 > 0.99);

	helper.chiSqFisher = chi_squared(helper[0],helper[1],helper[4],helper[5]);
	//helper.chiSqFisher = 0;

	return true;
}

template <class T>
bool LDSlave<T>::CalculateLDUnphasedMath(){
	// Total amount of non-missing alleles
	helper.totalHaplotypeCounts  = helper[0]  + helper[1]  + helper[4]  + helper[5]
	                          + helper[16] + helper[17] + helper[20] + helper[21]
	                          + helper[64] + helper[65] + helper[68] + helper[69]
	                          + helper[80] + helper[81] + helper[84] + helper[85];

	// All values are missing or too few
	if(helper.totalHaplotypeCounts < MINIMUM_ALLOWED_ALLELES){
		//++this->insufficent_alleles;
		return false;
	}


	// How many hets-hets is there? 0/1 0/1 or 1/0 1/0 or equivalent
	const float number_of_hets = helper[17] + helper[20] + helper[65] + helper[68];

	// If het-hets are 0 then do normal calculations
	// There is no phase uncertainty
	// Use phased math
	if(number_of_hets == 0){
		const float p0 = 2*helper[0] + helper[1]  + helper[4]    + helper[16] + helper[64];
		const float q0 = helper[16]  + helper[64] + 2*helper[80] + helper[81] + helper[84];
		const float p1 = helper[1]   + helper[4]  + 2*helper[5]  + helper[21] + helper[69];
		const float q1 = helper[21]  + helper[69] + helper[81]   + helper[84] + 2*helper[85];

		helper[0] = p0;
		helper[1] = p1;
		helper[4] = q0;
		helper[5] = q1;

		// Update counter
		//++this->no_uncertainty;

		// Reset
		helper.chiSqModel = 0;

		// Use standard math
		return(this->CalculateLDPhasedMath());
	}

	const double p = ((helper[0] + helper[1]  + helper[4]  + helper[5])*2.0
				   + (helper[16] + helper[17] + helper[20] + helper[21] + helper[64] + helper[65] + helper[68] + helper[69]))
				   / (2.0 * helper.totalHaplotypeCounts);
	const double q = ((helper[0] + helper[16] + helper[64] + helper[80])*2.0
				   + (helper[1]  + helper[4]  + helper[17] + helper[20] + helper[65] + helper[68] + helper[81] + helper[84]))
				   / (2.0 * helper.totalHaplotypeCounts);
	const double n11 = 2.0* helper[0] + helper[1] + helper[4] + helper[16] + helper[64];

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
		const double theta 		= acos(-yN / sqrt(h2)) / 3.0;
		const double constant 	= 2.0 * sqrt(d2);
		const double alpha 		= xN + constant * cos(theta);
		const double beta 		= xN + constant * cos(2.0*M_PI/3.0 + theta);
		const double gamma 		= xN + constant * cos(4.0*M_PI/3.0 + theta);

		BYTE biologically_possible = 0;
		helper.chiSqModel = std::numeric_limits<float>::max();
		const double* chosen = &alpha;
		if(alpha >= minhap - constants::ALLOWED_ROUNDING_ERROR && alpha <= maxhap + constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			helper.chiSqModel = this->ChiSquaredUnphasedTable(alpha, p, q);
			chosen = &alpha;
		}

		if(beta >= minhap - constants::ALLOWED_ROUNDING_ERROR && beta <= maxhap + constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			if(this->ChiSquaredUnphasedTable(beta, p, q) < helper.chiSqModel){
				chosen = &beta;
				helper.chiSqModel = this->ChiSquaredUnphasedTable(beta, p, q);
			}
		}

		if(gamma >= minhap - constants::ALLOWED_ROUNDING_ERROR && gamma <= maxhap + constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			if(this->ChiSquaredUnphasedTable(gamma, p, q) < helper.chiSqModel){
				chosen = &gamma;
				helper.chiSqModel = this->ChiSquaredUnphasedTable(gamma, p, q); // implicit
			}
		}

		if(biologically_possible == 0){
			//++this->impossible;
			return false;
		}

		if(biologically_possible > 1)
			helper.setMultipleRoots();

		//++this->possible;
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
		if(!(alpha >= minhap - constants::ALLOWED_ROUNDING_ERROR && alpha <= maxhap + constants::ALLOWED_ROUNDING_ERROR)){
			//++this->impossible;
			return false;
		}

		helper.chiSqModel = this->ChiSquaredUnphasedTable(alpha, p, q);
		//++this->possible;

		return(this->ChooseF11Calculate(alpha, p, q));

	} else { // Yn2 == h2
		const double delta = pow((yN/2.0*a),(1.0/3.0));
		const double alpha = xN + delta; // alpha = beta in this case
		const double gamma = xN - 2.0*delta;

		if(std::isnan(alpha) || std::isnan(gamma)){
			//++this->impossible;
			return false;
		}

		BYTE biologically_possible = 0;
		helper.chiSqModel = std::numeric_limits<float>::max();
		const double* chosen = &alpha;
		if(alpha >= minhap - constants::ALLOWED_ROUNDING_ERROR && alpha <= maxhap + constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			helper.chiSqModel = this->ChiSquaredUnphasedTable(alpha, p, q);
			chosen = &alpha;
		}

		if(gamma >= minhap - constants::ALLOWED_ROUNDING_ERROR && gamma <= maxhap + constants::ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			if(this->ChiSquaredUnphasedTable(gamma, p, q) < helper.chiSqModel){
				chosen = &gamma;
				helper.chiSqModel = this->ChiSquaredUnphasedTable(gamma, p, q); // implicit
			}
		}

		if(biologically_possible == 0){
			//++this->impossible;
			return false;
		}

		//++this->possible;
		return(this->ChooseF11Calculate(*chosen, p, q));
	}
	return(false);
}

template <class T>
bool LDSlave<T>::CalculateLDUnphasedVectorizedNoMissing(const block_type& block1, const block_type& block2){
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

	const simd_pair& datA = block1.currentBitvector();
	const simd_pair& datB = block2.currentBitvector();
	const BYTE* const arrayA = datA.data;
	const BYTE* const arrayB = datB.data;

#if SIMD_AVAILABLE == 1
	const U32 frontSmallest = datA.frontZero < datB.frontZero ? datA.frontZero : datB.frontZero;
	const U32 tailSmallest  = datA.tailZero  < datB.tailZero  ? datA.tailZero  : datB.tailZero;
	U32 i = frontSmallest;
	const U32 frontBonus    = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus     = (datA.tailZero != tailSmallest  ? datA.tailZero  : datB.tailZero);

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

	helper[0]  = helper_simd.counters[0] - this->unphased_unbalanced_adjustment;
	helper[0] += (frontSmallest + tailSmallest) * GENOTYPE_TRIP_COUNT;
	helper[1]  = helper_simd.counters[1];
	helper[5]  = helper_simd.counters[2];
	helper[16] = helper_simd.counters[3];
	helper[21] = helper_simd.counters[5];
	helper[80] = helper_simd.counters[6];
	helper[81] = helper_simd.counters[7];
	helper[85] = helper_simd.counters[8];
	//helper[17] = helper_simd.counters[4];
	helper[17] = this->samples - (helper[0] +  helper[1] + helper[5] + helper[16] + helper[21] + helper[80] + helper[81] + helper[85]);

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << ticks_per_iter.count() << '\n';
#endif

	this->setFLAGs(block1, block2);
	return(this->CalculateLDUnphasedMath());
}

template <class T>
bool LDSlave<T>::CalculateLDUnphasedVectorized(const block_type& block1, const block_type& block2){
	#if SLAVE_DEBUG_MODE < 6
	if(block1.currentMeta().has_missing == 0 && block2.currentMeta().has_missing == 0)
		return(this->CalculateLDUnphasedVectorizedNoMissing(block1, block2));
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

	const simd_pair& datA         = block1.currentBitvector();
	const simd_pair& datB         = block2.currentBitvector();
	const BYTE* const arrayA      = datA.data;
	const BYTE* const arrayB      = datB.data;
	const BYTE* const arrayA_mask = datA.mask;
	const BYTE* const arrayB_mask = datB.mask;

#if SIMD_AVAILABLE == 1
	const U32 frontSmallest = datA.frontZero < datB.frontZero ? datA.frontZero : datB.frontZero;
	const U32 tailSmallest  = datA.tailZero  < datB.tailZero  ? datA.tailZero  : datB.tailZero;
	U32 i = frontSmallest;
	const U32 frontBonus    = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus     = (datA.tailZero != tailSmallest  ? datA.tailZero  : datB.tailZero);

	//std::cerr << frontSmallest << '\t' << tailSmallest << std::endl;

	const VECTOR_TYPE* const vectorA      = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB      = (const VECTOR_TYPE* const)arrayB;
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
	POPCOUNT(helper_simd.counters[0], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);	\
	POPCOUNT(helper_simd.counters[1], __intermediate);				\
	__intermediate = FILTER_UNPHASED(altref, altref); 						\
	POPCOUNT(helper_simd.counters[2], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);	\
	POPCOUNT(helper_simd.counters[3], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refref, altalt, altalt, refref);	\
	POPCOUNT(helper_simd.counters[4], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altref, altref, refalt);	\
	POPCOUNT(helper_simd.counters[4], __intermediate);				\
	__intermediate = FILTER_UNPHASED(refalt,refalt);						\
	POPCOUNT(helper_simd.counters[6], __intermediate);				\
}

#define ITER_LONG {															\
	ITER_SHORT																\
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref); 	\
	POPCOUNT(helper_simd.counters[5], __intermediate);				\
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt); 	\
	POPCOUNT(helper_simd.counters[7], __intermediate);				\
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt); 						\
	POPCOUNT(helper_simd.counters[8], __intermediate);				\
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

	helper[0]  = helper_simd.counters[0] - this->unphased_unbalanced_adjustment;
	helper[0] += (frontSmallest + tailSmallest) * GENOTYPE_TRIP_COUNT;
	helper[1]  = helper_simd.counters[1];
	helper[5]  = helper_simd.counters[2];
	helper[16] = helper_simd.counters[3];
	helper[17] = helper_simd.counters[4];
	helper[21] = helper_simd.counters[5];
	helper[80] = helper_simd.counters[6];
	helper[81] = helper_simd.counters[7];
	helper[85] = helper_simd.counters[8];

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << ticks_per_iter.count() << '\n';
#endif

	this->setFLAGs(block1, block2);
	return(this->CalculateLDUnphasedMath());
}

template <class T>
bool LDSlave<T>::CalculateLDPhasedVectorized(const block_type& block1, const block_type& block2){
	#if SLAVE_DEBUG_MODE < 6
	if(block1.currentMeta().has_missing == 0 && block2.currentMeta().has_missing == 0){
		return(this->CalculateLDPhasedVectorizedNoMissing(block1, block2));
		//return(this->CalculateLDPhasedVectorizedNoMissingNoTable(helper, block1, block2));
	}
	#endif

	helper.resetPhased();
	helper_simd.counters[0] = 0;
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;

	const simd_pair& datA         = block1.currentBitvector();
	const simd_pair& datB         = block2.currentBitvector();
	const BYTE* const arrayA      = datA.data;
	const BYTE* const arrayB      = datB.data;
	const BYTE* const arrayA_mask = datA.mask;
	const BYTE* const arrayB_mask = datB.mask;

#if SIMD_AVAILABLE == 1
	const U32 frontSmallest = datA.frontZero < datB.frontZero ? datA.frontZero : datB.frontZero;
	const U32 tailSmallest  = datA.tailZero  < datB.tailZero  ? datA.tailZero  : datB.tailZero;
	U32 i = frontSmallest;
	const U32 frontBonus    = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus     = (datA.tailZero != tailSmallest  ? datA.tailZero  : datB.tailZero);

	const VECTOR_TYPE* const vectorA      = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB      = (const VECTOR_TYPE* const)arrayB;
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
	POPCOUNT(helper_simd.counters[0], __intermediate);				\
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(helper_simd.counters[1], __intermediate);				\
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(helper_simd.counters[2], __intermediate);				\
	i += 1;																	\
}

#define ITER {																\
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);					\
	__intermediate  = PHASED_ALTALT_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(helper_simd.counters[3], __intermediate);				\
	__intermediate  = PHASED_REFREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(helper_simd.counters[0], __intermediate);				\
	__intermediate  = PHASED_ALTREF_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(helper_simd.counters[1], __intermediate);				\
	__intermediate  = PHASED_REFALT_MASK(vectorA[i], vectorB[i], masks);	\
	POPCOUNT(helper_simd.counters[2], __intermediate);				\
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
			helper_simd.scalarA[l] = (arrayA[k+l] & arrayB[k+l]) & mask;
			helper_simd.scalarB[l] = ((~arrayA[k+l]) & (~arrayB[k+l])) & mask;
			helper_simd.scalarC[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayA[k+l]) & mask;
			helper_simd.scalarD[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayB[k+l]) & mask;
		}
		helper_simd.counters[0] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(helper_simd.scalarB));
		helper_simd.counters[2] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(helper_simd.scalarC));
		helper_simd.counters[1] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(helper_simd.scalarD));
		helper_simd.counters[3] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(helper_simd.scalarA));
	}

#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k < this->byte_width; ++k){
		mask = ~(arrayA_mask[k] | arrayB_mask[k]);
		helper_simd.counters[0] += POPCOUNT_ITER(((~arrayA[k]) & (~arrayB[k])) & mask);
		helper_simd.counters[2] += POPCOUNT_ITER(((arrayA[k] ^ arrayB[k]) & arrayA[k]) & mask);
		helper_simd.counters[1] += POPCOUNT_ITER(((arrayA[k] ^ arrayB[k]) & arrayB[k]) & mask);
		helper_simd.counters[3] += POPCOUNT_ITER((arrayA[k] & arrayB[k]) & mask);
	}

	helper[1] = helper_simd.counters[1];
	helper[4] = helper_simd.counters[2];
	helper[5] = helper_simd.counters[3];
	helper[0] = (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 + helper_simd.counters[0] - this->phased_unbalanced_adjustment;


#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << "V\t" << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif

	this->setFLAGs(block1, block2);
	return(this->CalculateLDPhasedMath());
}

template <class T>
bool LDSlave<T>::CalculateLDPhasedVectorizedNoMissing(const block_type& block1, const block_type& block2){
	helper.resetPhased();
	helper_simd.counters[0] = 0;
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;

	const simd_pair& datA    = block1.currentBitvector();
	const simd_pair& datB    = block2.currentBitvector();
	const BYTE* const arrayA = datA.data;
	const BYTE* const arrayB = datB.data;

#if SIMD_AVAILABLE == 1
	const U32 frontSmallest = datA.frontZero < datB.frontZero ? datA.frontZero : datB.frontZero;
	const U32 tailSmallest  = datA.tailZero  < datB.tailZero  ? datA.tailZero  : datB.tailZero;
	U32 i = frontSmallest;
	const U32 frontBonus    = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus     = datA.tailZero  != tailSmallest  ? datA.tailZero  : datB.tailZero;

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
	POPCOUNT(helper_simd.counters[2], __intermediate);	\
	__intermediate  = PHASED_ALTREF(vectorA[i], vectorB[i]);	\
	POPCOUNT(helper_simd.counters[1], __intermediate);	\
	i += 1;														\
}

#define ITER {													\
	__intermediate  = PHASED_ALTALT(vectorA[i], vectorB[i]);	\
	POPCOUNT(helper_simd.counters[3], __intermediate);	\
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
			helper_simd.scalarA[l] = (arrayA[k+l] & arrayB[k+l]);
			helper_simd.scalarC[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayA[k+l]);
			helper_simd.scalarD[l] = ((arrayA[k+l] ^ arrayB[k+l]) & arrayB[k+l]);
		}

		helper_simd.counters[2] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(helper_simd.scalarC));
		helper_simd.counters[1] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(helper_simd.scalarD));
		helper_simd.counters[3] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(helper_simd.scalarA));
	}

#ifdef __INTEL_COMPILER
	#pragma vector aligned
#endif
	for(; k < this->byte_width; ++k){
		helper_simd.counters[2] += POPCOUNT_ITER((arrayA[k] ^ arrayB[k]) & arrayA[k]);
		helper_simd.counters[1] += POPCOUNT_ITER((arrayA[k] ^ arrayB[k]) & arrayB[k]);
		helper_simd.counters[3] += POPCOUNT_ITER(arrayA[k] & arrayB[k]);
	}

	helper[1] = helper_simd.counters[1];
	helper[4] = helper_simd.counters[2];
	helper[5] = helper_simd.counters[3];
	helper[0] = this->samples*2 - (helper[1] + helper[4] + helper[5]);
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << "V\t" << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << "\t" << ticks_per_iter.count() << '\n';
#endif

	this->setFLAGs(block1, block2);

	return(this->CalculateLDPhasedMath());
}

template <class T>
bool LDSlave<T>::CalculateLDPhasedVectorizedNoMissingNoTable(const block_type& block1, const block_type& block2){
	helper.resetPhased();
	helper_simd.counters[0] = 0;

	const simd_pair& datA    = block1.currentBitvector();
	const simd_pair& datB    = block2.currentBitvector();
	const BYTE* const arrayA = datA.data;
	const BYTE* const arrayB = datB.data;

#if SIMD_AVAILABLE == 1
	const U32 frontSmallest = datA.frontZero < datB.frontZero ? datA.frontZero : datB.frontZero;
	const U32 tailSmallest  = datA.tailZero  < datB.tailZero  ? datA.tailZero  : datB.tailZero;
	const U32 frontBonus    = datA.frontZero != frontSmallest ? datA.frontZero : datB.frontZero;
	const U32 tailBonus     = datA.tailZero  != tailSmallest  ? datA.tailZero  : datB.tailZero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	VECTOR_TYPE __intermediate;

#define ITER_SHORT {											\
	__intermediate  = PHASED_REFREF(vectorA[i], vectorB[i]);	\
	POPCOUNT(helper_simd.counters[0], __intermediate);			\
	i += 1;														\
}

	U32 i = frontSmallest;
	for( ; i < frontBonus; ) 					  	ITER_SHORT
	for( ; i < this->vectorCycles - tailBonus; )  	ITER_SHORT
	for( ; i  < this->vectorCycles - tailSmallest; )ITER_SHORT

#undef ITER_SHORT
	U32 k = this->byteAlignedEnd;
#else
	U32 k = 0;
#endif

	for(; k+8 < this->byte_width; k += 8){
		for(U32 l = 0; l < 8; ++l)
			helper_simd.scalarB[l] = (~arrayA[k+l]) & (~arrayB[k+l]);
		helper_simd.counters[0] += POPCOUNT_ITER(*reinterpret_cast<const U64* const>(helper_simd.scalarB));
	}

	for(; k < this->byte_width; ++k){
		helper_simd.counters[0] += POPCOUNT_ITER(((~arrayA[k]) & (~arrayB[k]) & 255));
	}
	helper[0] = (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 + helper_simd.counters[0] - this->phased_unbalanced_adjustment;

	//this->setFLAGs(block1, block2);
	return(this->CalculateLDPhasedMathSimple(block1, block2));
}

template <class T>
bool LDSlave<T>::CalculateLDPhasedSimple(const block_type& block1, const block_type& block2){
	helper.resetPhased();


	const base::HaplotypeBitVector& bitvectorA = block1.currentHaplotypeBitvector();
	const base::HaplotypeBitVector& bitvectorB = block2.currentHaplotypeBitvector();

	const U32 n_cycles = std::min(bitvectorA.l_list, bitvectorB.l_list);
	const U32 n_total = bitvectorA.l_list + bitvectorB.l_list;
	U32 n_same = 0;

	// compare A to B
	if(bitvectorA.l_list >= bitvectorB.l_list){
		for(U32 i = 0; i < n_cycles; ++i){
			n_same += bitvectorA.get(bitvectorB.indices[i]);
		}
	} else {
		for(U32 i = 0; i < n_cycles; ++i){
			n_same += bitvectorB.get(bitvectorA.indices[i]);
		}
	}

	helper[0] = 2*this->samples - (n_total - n_same);

	return(this->CalculateLDPhasedMathSimple(block1, block2));
}

template <class T>
bool LDSlave<T>::CalculateLDPhasedSimpleSample(const block_type& block1, const block_type& block2){
	helper.resetPhased();

	const base::HaplotypeBitVector& bitvectorA = block1.currentHaplotypeBitvector();
	const base::HaplotypeBitVector& bitvectorB = block2.currentHaplotypeBitvector();

	const U32 n_cycles = 1000;
	const int64_t n_total = bitvectorA.l_list + bitvectorB.l_list; // Number of positions that are non-reference
	int64_t n_same = 0;

	// compare A to B
	if(bitvectorA.l_list >= bitvectorB.l_list){
		for(U32 i = 0; i < n_cycles; ++i){
			n_same += bitvectorA.get(bitvectorB.sampled_indicies[i]); // If bit[position] is set in both then is not jointly REF
		}
		// Linear estimator
		n_same *= (double)bitvectorB.l_list / n_cycles;
	} else {
		for(U32 i = 0; i < n_cycles; ++i){
			n_same += bitvectorB.get(bitvectorA.sampled_indicies[i]); // If bit[position] is set in both then is not jointly REF
		}
		// Linear estimator
		n_same *= (double)bitvectorA.l_list / n_cycles;
	}

	// Lower bounds
	n_same = std::max(n_total - 2*(int64_t)this->samples, n_same);

	/*
	assert(n_same <= n_total);
	if(2*(int64_t)this->samples - (n_total - n_same) < 0){
		std::cerr << 2*(int64_t)this->samples - (n_total - n_same) << "\t" << n_total << "\t" << n_same << "\t" << n_total - n_same << std::endl;
		std::cerr << n_total - 2*this->samples << " OR " << n_same << std::endl;
		std::cerr << std::max(n_total - 2*this->samples, n_same) << std::endl;
		exit(1);
	}
	*/

	//std::cerr << n_same << "\t" << factor << "\t" << n_same*factor << "\t" << block1.currentMeta().AF << "\t" << block2.currentMeta().AF << std::endl;

	// Number of joint REF = 2N - (total alleles not REF - total joint REF)
	helper[0] = 2*this->samples - (n_total - n_same);
	helper.setSampled(true);

	return(this->CalculateLDPhasedMathSimple(block1, block2));
	/*
	{
		std::cerr << helper[0] << ": " << helper.R2 << " with " << n_same << std::endl;
	}
	*/
}

template <class T>
bool LDSlave<T>::CalculateLDPhased(const block_type& block1, const block_type& block2){
	helper.resetPhased();
#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE == 5
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	const run_type* const a = block1.current();
	const run_type* const b = block2.current();
	T currentLengthA = a[0].runs;
	T currentLengthB = b[0].runs;

	BYTE currentMixL = ((  a[0].alleleA & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH)
						^ (b[0].alleleA & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1));
	BYTE currentMixR = ((  a[0].alleleB & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH)
						^ (b[0].alleleB & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1));

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
		} else if(currentLengthA < currentLengthB){ // If processed run length A < processed run length B
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
		helper[currentMixL] += add;
		helper[currentMixR] += add;

		// Exit condition
		if(pointerA == block1.currentMeta().runs || pointerB == block2.currentMeta().runs){
			if(pointerA != block1.currentMeta().runs || pointerB != block2.currentMeta().runs){
				std::cerr << helpers::timestamp("FATAL") << "Failed to exit equally!\n" << pointerA << "/" << block1.currentMeta().runs << " and " << pointerB << "/" << block2.currentMeta().runs << std::endl;
				exit(1);
			}
			break;
		}

		// Update mixing value
		currentMixL = ((  a[pointerA].alleleA & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH)
					   ^ (b[pointerB].alleleA & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1));
		currentMixR = ((  a[pointerA].alleleB & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1)) << constants::TOMAHAWK_ALLELE_PACK_WIDTH)
					   ^ (b[pointerB].alleleB & ((1 << constants::TOMAHAWK_ALLELE_PACK_WIDTH)-1));

	#if SLAVE_DEBUG_MODE == 3
		++iterations;
		std::cout << a.getMeta().runs << '\t' << b.getMeta().runs << '\t' << iterations << std::endl;
	#endif
	}

#if SLAVE_DEBUG_MODE == 4 || SLAVE_DEBUG_MODE ==  5
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	std::cout << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif

	this->setFLAGs(block1, block2);

	return(this->CalculateLDPhasedMath());
}

template <class T>
bool LDSlave<T>::CalculateLDPhasedMath(){
	// Trigger phased flag
	helper.setUsedPhasedMath();

	// Total amount of non-missing alleles
	helper.totalHaplotypeCounts = helper[0] + helper[1] + helper[4] + helper[5];

	// All values are missing
	if(helper.totalHaplotypeCounts < MINIMUM_ALLOWED_ALLELES){
		//++this->insufficent_alleles;
		return false;
	}
	//++this->possible;

	// Filter by total minor haplotype frequency
	if(helper.getTotalAltHaplotypeCount() < this->parameters.minimum_sum_alternative_haplotype_count){
		//std::cerr << "FILTER: " << helper.getTotalAltHaplotypeCount() << std::endl;
		return false;
	}


	// Haplotype frequencies
	helper.haplotypeCounts[0] = (helper[0] + helper[1]) / helper.totalHaplotypeCounts; // pA
	helper.haplotypeCounts[1] = (helper[4] + helper[5]) / helper.totalHaplotypeCounts; // qA
	helper.haplotypeCounts[2] = (helper[0] + helper[4]) / helper.totalHaplotypeCounts; // pB
	helper.haplotypeCounts[3] = (helper[1] + helper[5]) / helper.totalHaplotypeCounts; // qB

	const double divisor = helper.haplotypeCounts[0]*helper.haplotypeCounts[1]*helper.haplotypeCounts[2]*helper.haplotypeCounts[3];
	if(divisor == 0) return false;

	helper.D  = helper[0]/helper.totalHaplotypeCounts * helper[5]/helper.totalHaplotypeCounts - helper[1]/helper.totalHaplotypeCounts * helper[4]/helper.totalHaplotypeCounts;
	helper.R2 = helper.D*helper.D / divisor;
	helper.R  = sqrt(helper.R2);

	if(helper.R2 >= this->parameters.R2_min && helper.R2 <= this->parameters.R2_max){
		if(helper.D >= 0){
			helper.Dmax = helper.haplotypeCounts[0]*helper.haplotypeCounts[3] < helper.haplotypeCounts[1]*helper.haplotypeCounts[2]
					? helper.haplotypeCounts[0]*helper.haplotypeCounts[3]
					: helper.haplotypeCounts[1]*helper.haplotypeCounts[2];
		} else {
			helper.Dmax = helper.haplotypeCounts[0]*helper.haplotypeCounts[2] < helper.haplotypeCounts[1]*helper.haplotypeCounts[3]
					? -helper.haplotypeCounts[0]*helper.haplotypeCounts[2]
					: -helper.haplotypeCounts[1]*helper.haplotypeCounts[3];
		}
		helper.Dprime = helper.D / helper.Dmax;

		// Calculate Fisher's exact test P-value
		double left,right,both;
		kt_fisher_exact(round(helper[0]),round(helper[1]),
		                round(helper[4]),round(helper[5]),
						&left, &right, &both);
		helper.P = both;

		// Fisher's exact test P value filter
		if(helper.P > this->parameters.P_threshold){
			return false;
		}

		helper.setCompleteLD(helper[0] == 0 || helper[1] == 0 || helper[4] == 0 || helper[5] == 0);
		helper.setPerfectLD(helper.R2 > 0.99);

		// Calculate Chi-Sq CV from 2x2 contingency table
		helper.chiSqModel = 0;
		//helper.chiSqFisher = chi_squared(helper[0],helper[1],helper[4],helper[5]);
		helper.chiSqFisher = helper.totalHaplotypeCounts * helper.R2;
		//helper.chiSqFisher = 0;

		return true;
	}
	return false;
}

template <class T>
bool LDSlave<T>::CalculateLDPhasedMathSimple(const block_type& block1, const block_type& block2){
	//++this->possible;

	// If haplotype count for (0,0) >
	if(helper[0] > 2*this->samples - this->parameters.minimum_sum_alternative_haplotype_count)
		return false;

	// D = (joint HOM_HOM) - (HOM_A * HOM_B) = pAB - pApB
	const double divisor = block1.currentMeta().AF * (1 - block1.currentMeta().AF) * block2.currentMeta().AF * (1 - block2.currentMeta().AF);
	if(divisor == 0) return false;

	helper.D = helper[0]/(2*this->samples) - block1.currentMeta().AF*block2.currentMeta().AF;
	helper.R2 = helper.D*helper.D / divisor ;

	if(helper.R2 >= this->parameters.R2_min && helper.R2 <= this->parameters.R2_max){
		helper.R  = sqrt(helper.R2);

		if(helper.D < 0){
			if( block1.currentMeta().AF * block2.currentMeta().AF < (1 - block1.currentMeta().AF) * (1 - block2.currentMeta().AF)){
				helper.Dprime = helper.D / (block1.currentMeta().AF * block2.currentMeta().AF);
			} else {
				helper.Dprime = helper.D / ( (1 - block1.currentMeta().AF) * (1 - block2.currentMeta().AF) );
			}
		} else { // D >= 0
			if((1-block1.currentMeta().AF) * block2.currentMeta().AF < block1.currentMeta().AF * (1 - block2.currentMeta().AF)){
				helper.Dprime = helper.D / ((1 - block1.currentMeta().AF) * block2.currentMeta().AF);
			} else {
				helper.Dprime = helper.D / (block1.currentMeta().AF * (1 - block2.currentMeta().AF));
			}
		}

		// Trigger phased flag
		helper.setUsedPhasedMath();
		helper.setFastMode();
		helper.setPerfectLD(helper.R2 > 0.99);
		helper.setCompleteLD(helper.Dprime > 0.99);
		this->setFLAGs(block1, block2);

		return true;
	}

	return false;
}

// Execute diagonal working order
template <class T>
bool LDSlave<T>::DiagonalWorkOrder(const order_type& order){
	block_type block1(this->manager[order.row]);

	for(U32 j = order.fromColumn; j < order.toColumn; ++j){
		//std::cerr << helpers::timestamp("DEBUG", "DIAG") << i << '/' << j << '\t' << order << std::endl;
		if(order.row == j)
			this->CompareBlocks(helper, block1);
		else {
			block_type block2(this->manager[j]);
			this->CompareBlocks(block1, block2);
		}
	}

	return true;
}

// Execute square working order
template <class T>
bool LDSlave<T>::SquareWorkOrder(const order_type& order){
	block_type block1(this->manager[order.row]);

	for(U32 j = order.fromColumn; j < order.toColumn; ++j){
		//std::cerr << helpers::timestamp("DEBUG", "SQUARE") << i << '/' << j << '\t' << order << std::endl;
		if(order.row == j)
			this->CompareBlocks(block1);
		else {
			block_type block2(this->manager[j]);
			this->CompareBlocks(block1, block2);
		}
	}

	return true;
}

template <class T>
bool LDSlave<T>::Calculate(void){
	// If there is no data
	if(this->manager.size() == 0){
		std::cerr << helpers::timestamp("ERROR", "CONTROLLER") << "There is no data..." << std::endl;
		return false;
	}

	// If there is no data
	if(this->orders.size() == 0){
		std::cerr << helpers::timestamp("WARNING", "SLAVE") << "No data was given to this worker. Balancing incomplete..." << std::endl;
		return false;
	}

	// Foreach work order
	for(U32 i = 0; i < this->orders.size(); ++i){
		// If the order includes part of the diagonal
		if(!this->SquareWorkOrder(this->orders[i]))
			return false;
	}

	// Finish the output manager
	this->output_writer.flush();

	return true;
}

template <class T>
bool LDSlave<T>::CompareBlocksFunction(const block_type& block1, const block_type& block2){
#if SLAVE_DEBUG_MODE == 1 // 1 = No debug mode
	if(block1.currentMeta().all_phased == 1 && block2.currentMeta().all_phased == 1){
		if(block1.currentMeta().runs + block2.currentMeta().runs <= 20){
			return(this->CalculateLDPhased(block1, block2));
		} else {
			return(this->CalculateLDPhasedVectorized(block1, block2));
		}
	} else {
		if(block1.currentMeta().runs + block2.currentMeta().runs <= 20){
			return(this->CalculateLDUnphased(block1, block2));
		} else {
			return(this->CalculateLDUnphasedVectorized(block1, block2));
		}
	}

#elif SLAVE_DEBUG_MODE == 2
	if(this->CalculateLDPhasedVectorizedNoMissing(block1, block2)){
		this->output_writer.Add(block1, block2, helper);
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
	m.addPhased(0, helper);
	this->CalculateLDPhasedVectorized(block1, block2);
	m.addPhased(1, helper);
	this->CalculateLDPhasedVectorizedNoMissing(block1, block2);
	m.addPhased(2, helper);
	this->CalculateLDUnphased(block1, block2);
	m.addUnphased(0, helper);
	this->CalculateLDUnphasedVectorized(block1, block2);
	m.addUnphased(1, helper);
	this->CalculateLDUnphasedVectorizedNoMissing(block1, block2);
	m.addUnphased(2, helper);

	if(!m.validate()){
		std::cerr << helpers::timestamp("ERROR", "VALIDATION") << "Failed validation..." << std::endl;
		std::cerr << m << std::endl;
		//exit(1);
	}
#elif SLAVE_DEBUG_MODE == 8
	this->CalculateLDPhasedVectorizedNoMissingNoTable(block1, block2);
	const double r1 = helper.R2;
	const double c1 = helper[0];

	this->CalculateLDPhasedVectorizedNoMissing(block1, block2);
	const double r2 = helper.R2;
	const double c2 = helper[0];

	if(abs(r1-r2) > 0.05){
		std::cerr << r1 << "\t" << r2 << "\t" << block1.currentMeta().AF << "\t" << block2.currentMeta().AF << '\t' << helper[0] << "," << helper[1] << "," << helper[4] << "," << helper[5] << std::endl;
	}
#endif
}

template <class T>
bool LDSlave<T>::CompareBlocksFunctionForcedPhased(const block_type& block1, const block_type& block2){
	if(block1.currentMeta().runs + block2.currentMeta().runs <= 200){
		return(this->CalculateLDPhased(block1, block2));
	} else {
		return(this->CalculateLDPhasedVectorized(block1, block2));
	}
}

template <class T>
bool LDSlave<T>::CompareBlocksFunctionForcedPhasedSimple(const block_type& block1, const block_type& block2){
	if(std::min(block1.currentHaplotypeBitvector(). l_list,block2.currentHaplotypeBitvector().l_list) <= 1000){
		if(block1.currentMeta().has_missing == false && block2.currentMeta().has_missing == false)
			return(this->CalculateLDPhasedSimple(block1, block2));
		return(this->CalculateLDPhased(block1, block2));
	} else {
		if(block1.currentMeta().has_missing == false && block2.currentMeta().has_missing == false)
			return(this->CalculateLDPhasedVectorizedNoMissingNoTable(block1, block2));
		return(this->CalculateLDPhasedVectorized(block1, block2));
	}
}

template <class T>
bool LDSlave<T>::CompareBlocksFunctionForcedUnphased(const block_type& block1, const block_type& block2){
	// Ignore when one or both is invariant
	if(block1.currentMeta().AF  == 0 || block2.currentMeta().AF  == 0 ||
       block1.currentMeta().runs == 1 || block2.currentMeta().runs == 1)
	{
		//std::cerr << "invariant" << std::endl;
		return false;
	}

	if(block1.currentMeta().runs + block2.currentMeta().runs <= 200){
		return(this->CalculateLDUnphased(block1, block2));
	} else {
		return(this->CalculateLDUnphasedVectorized(block1, block2));
	}
}

// Within-block comparisons
template <class T>
bool LDSlave<T>::CompareBlocks(block_type& block1){
	//std::cerr << helpers::timestamp("DEBUG", "DIAG-INTERNAL") << (block1.size()*block1.size()-block1.size())/2 << std::endl;
	block1.resetIterator(); // make sure it is reset
	block_type block2(block1);

	for(U32 i = 0; i < block1.size(); ++i){
 		block2 = block1;
		++block2; // block2 starts at relative +1

		for(U32 j = i + 1; j < block1.size(); ++j){
			// temp
			//if(block2.currentMeta().position - block1.currentMeta().position > 1000000)
			//	continue;

			if((this->*phase_function_across)(block1, block2)){
				this->output_writer.Add(block1.currentMeta(), block2.currentMeta(), block1.getTotempole(), block2.getTotempole(), helper);
			}
			++block2;
		}

		// Update progress
		this->progress(block1.size() - (i + 1), this->output_writer.getProgressCounts());
		this->output_writer.ResetProgress();
		++block1;
	}
	return true;
}

// Across block comparisons
template <class T>
bool LDSlave<T>::CompareBlocks(block_type& block1, block_type& block2){
	//std::cerr << helpers::timestamp("DEBUG", "DIAG-SQUARE") << block1.size()*block2.size() << std::endl;

	// Reset
	block1.resetIterator();
	block2.resetIterator();

	//if(block2.getMeta(0).position - block1.getMeta(0).position > 1000000)
	//	return false;

	// Cycle over block 1 and block 2
	for(U32 i = 0; i < block1.size(); ++i){
		for(U32 j = 0; j < block2.size(); ++j){
			// temp
			//if(block2.currentMeta().position - block1.currentMeta().position > 1000000)
			//	break;

			if((this->*phase_function_across)(block1, block2)){
				this->output_writer.Add(block1.currentMeta(), block2.currentMeta(), block1.getTotempole(), block2.getTotempole(), helper);
			}
			++block2;
		}

		// Update progress
		this->progress(block2.size(), this->output_writer.getProgressCounts());
		this->output_writer.ResetProgress();

		// Reset position in block2 and increment position in block1
		block2.resetIterator();
		++block1;
	}
	return true;
}

}

#endif /* TOMAHAWK_LD_CALCULATION_SLAVE_H_ */
