#ifndef ALGORITHM_TOMAHAWKSLAVESIMDHELPER_H_
#define ALGORITHM_TOMAHAWKSLAVESIMDHELPER_H_

#include "support/simd_definitions.h"

namespace tomahawk{
namespace support{

template <int Y = SIMD_ALIGNMENT>
struct LDCalculationSIMDHelper{
public:
	LDCalculationSIMDHelper(void) :
#if SIMD_AVAILABLE == 1
		counters((U64*)_mm_malloc(sizeof(U64)*16, Y)),
		scalarA((BYTE*)_mm_malloc(sizeof(BYTE)*8, Y)),
		scalarB((BYTE*)_mm_malloc(sizeof(BYTE)*8, Y)),
		scalarC((BYTE*)_mm_malloc(sizeof(BYTE)*8, Y)),
		scalarD((BYTE*)_mm_malloc(sizeof(BYTE)*8, Y))
#else
		counters(new U64[16]),
		scalarA(new BYTE[8]),
		scalarB(new BYTE[8]),
		scalarC(new BYTE[8]),
		scalarD(new BYTE[8])
#endif
	{
		memset(this->counters, 0, sizeof(U64)*16);
	}

	~LDCalculationSIMDHelper(){
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

public:
	U64*  counters;
	BYTE* scalarA;
	BYTE* scalarB;
	BYTE* scalarC;
	BYTE* scalarD;
} __attribute__((aligned(16)));

}
}

#endif /* ALGORITHM_TOMAHAWKSLAVESIMDHELPER_H_ */
