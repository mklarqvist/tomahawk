#ifndef TYPEDEFINITIONS_H_
#define TYPEDEFINITIONS_H_

//**************************************
// Basic Types
//**************************************
#if defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L   // C99
#include <stdint.h>
typedef int8_t		SBYTE;
typedef uint8_t		BYTE;
typedef int16_t		S16
typedef uint16_t	U16;
typedef uint32_t	U32;
typedef int32_t		S32;
typedef uint64_t	U64;
typedef uint64_t	ULL;
#else
typedef char				SBYTE;
typedef unsigned char		BYTE;
typedef unsigned short		U16;
typedef short				S16;
typedef unsigned int		U32;
typedef signed int			S32;
typedef unsigned long long	U64;
typedef unsigned long long	ULL;
#endif

#endif /* TYPEDEFINITIONS_H_ */
