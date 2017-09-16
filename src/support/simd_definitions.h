#ifndef SIMD_H_
#define SIMD_H_

#if defined(_MSC_VER)
     /* Microsoft C/C++-compatible compiler */
     #include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
     /* GCC-compatible compiler, targeting x86/x86-64 */
     #include <x86intrin.h>
#elif defined(__GNUC__) && defined(__ARM_NEON__)
     /* GCC-compatible compiler, targeting ARM with NEON */
     #include <arm_neon.h>
#elif defined(__GNUC__) && defined(__IWMMXT__)
     /* GCC-compatible compiler, targeting ARM with WMMX */
     #include <mmintrin.h>
#elif (defined(__GNUC__) || defined(__xlC__)) && (defined(__VEC__) || defined(__ALTIVEC__))
     /* XLC or GCC-compatible compiler, targeting PowerPC with VMX/VSX */
     #include <altivec.h>
#elif defined(__GNUC__) && defined(__SPE__)
     /* GCC-compatible compiler, targeting PowerPC with SPE */
     #include <spe.h>
#endif

const std::vector<std::string> SIMD_MAPPING = {"NONE","SSE","SSE2","SSE3","AVX","AVX2-256","AVX-512"};

#if defined(__AVX512F__) && __AVX512F__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	6
#define SIMD_ALIGNMENT	64
#define GENOTYPE_TRIP_COUNT	256
#elif defined(__AVX2__) && __AVX2__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	5
#define SIMD_ALIGNMENT	32
#define GENOTYPE_TRIP_COUNT	128
#elif defined(__AVX__) && __AVX__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	4
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#elif defined(__SSE4_1__) && __SSE4_1__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	3
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#elif defined(__SSE2__) && __SSE2__ == 1
#define SIMD_AVAILABLE	1
#define SIMD_VERSION	2
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#elif defined(__SSE__) && __SSE__ == 1
#define SIMD_AVAILABLE	0 // unsupported version
#define SIMD_VERSION	1
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#else
#define SIMD_AVAILABLE	0
#define SIMD_VERSION	0
#define SIMD_ALIGNMENT	16
#define GENOTYPE_TRIP_COUNT	64
#endif

#endif /* SIMD_H_ */
