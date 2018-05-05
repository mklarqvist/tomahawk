#ifndef TOMAHAWK_BASE_GENOTYPE_BITVECTOR_H_
#define TOMAHAWK_BASE_GENOTYPE_BITVECTOR_H_

#include "../tomahawk/genotype_objects.h"
#include "../support/simd_definitions.h"


namespace Tomahawk{
namespace Base{

/**<
 * Data structure to representing a 1-bit allele
 * representation of genotypes. This data structure
 * has to be aligned to `SIMD_ALIGNMENT` as specified
 * by the user CPU architecture. If no SIMD is available
 * on the device then use regular memory alignment.
 *
 * Special techniques to accelerate pairwise comparisons:
 * 1) Front and tail number of SIMD _elements_ (e.g. 128 bits / 16 bytes)
 *    that are either all 0 or 1. This allows the algorithm
 *    to either completely skip these stretches or
 *    resort to cheaper comparison functors.
 * 2) Counts of missingness needs to be maintained for these
 *    tail and head elements to function correctly.
 */
template <int T = SIMD_ALIGNMENT>
struct GenotypeBitvector{
public:
	GenotypeBitvector(const U32 size):
		frontZero(0),
		tailZero(0),
		frontZeroMissing(0),
		tailZeroMissing(0),
	#if SIMD_AVAILABLE == 1
		data((BYTE*)_mm_malloc(size, T)),
		mask((BYTE*)_mm_malloc(size, T))
	#else
		data(new BYTE[size]),
		mask(new BYTE[size])
	#endif
	{
		memset(this->data, 0, size);
		memset(this->mask, 0, size);
	}

	~GenotypeBitvector(){
	#if SIMD_AVAILABLE == 1
		_mm_free(this->data);
		_mm_free(this->mask);
	#else
		delete [] this->data;
		delete [] this->mask;
	#endif
	}

	template <class RLE_TYPE>
	bool build(const MetaEntry<RLE_TYPE>& meta_entry, const RLE_TYPE* const genotypes, const U32& n_samples){
		const U32 byte_width = ceil((double)n_samples/4);

		// INVERSE mask is cheaper in terms of instructions used
		// exploited in calculations: TomahawkCalculationSlave
		const BYTE lookup_mask[16] = {0, 0, 3, 3, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
		const BYTE lookup_data[16] = {0, 1, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

		Algorithm::GenotypeBitPacker packerA(this->data, 2);
		Algorithm::GenotypeBitPacker packerB(this->mask, 2);

		// Cycle over runs in container
		for(U32 j = 0; j < meta_entry.runs; ++j){
			const Support::GenotypeDiploidRunPacked<RLE_TYPE>* const packed = reinterpret_cast<const Support::GenotypeDiploidRunPacked<RLE_TYPE>* const>(&genotypes);
			packerA.add(lookup_data[packed->alleles], packed->runs);
			packerB.add(lookup_mask[packed->alleles], packed->runs);
		}

		const U32 byteAlignedEnd  = byte_width / (GENOTYPE_TRIP_COUNT/4) * (GENOTYPE_TRIP_COUNT/4);

		S32 j = 0;
		// Search from left->right
		for(; j < byteAlignedEnd; ++j){
			if(this->data[j] != 0 || this->mask[j] != 0)
				break;
		}

		// Front of zeroes
		this->frontZero = ((j - 1 < 0 ? 0 : j - 1)*4)/GENOTYPE_TRIP_COUNT;
		if(j != byteAlignedEnd){
			j = byteAlignedEnd - 1;
			for(; j > 0; --j){
				if(this->data[j] != 0 || this->mask[j] != 0)
					break;
			}
		}

		// Tail of zeroes
		this->tailZero = ((byteAlignedEnd - (j+1))*4)/GENOTYPE_TRIP_COUNT;

		return(true);
	}

public:
	U16   frontZero;        // leading zeros in aligned vector width
	U16   tailZero;         // trailing zeros in aligned vector width
	U16   frontZeroMissing; // number of missing values in leading zeros
	U16   tailZeroMissing;  // number of missing values in trailing zeros
	BYTE* data;
	BYTE* mask;
} __attribute__((aligned(16)));

}
}

#endif /* TOMAHAWK_BASE_GENOTYPE_BITVECTOR_H_ */
