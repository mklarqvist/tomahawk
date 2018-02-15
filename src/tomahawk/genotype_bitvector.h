#ifndef TOMAHAWK_BASE_GENOTYPE_BITVECTOR_H_
#define TOMAHAWK_BASE_GENOTYPE_BITVECTOR_H_

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
