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

}
