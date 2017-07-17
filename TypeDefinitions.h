/*
POLYGON - Next generation sequencing quality control and preprocessing
Copyright (C) 2015, Marcus D. R. Klarqvist.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#ifndef TYPEDEFINITIONS_H_
#define TYPEDEFINITIONS_H_

//**************************************
// Basic Types
//**************************************
#if defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L   // C99
#include <stdint.h>
typedef int8_t	SBYTE;
typedef uint8_t	BYTE;
typedef uint16_t	U16;
typedef uint32_t	U32;
typedef int32_t	S32;
typedef uint64_t	U64;
typedef uint64_t	ULL;
#else
typedef char 	SBYTE;
typedef unsigned char 	BYTE;
typedef unsigned short	U16;
typedef unsigned int		U32;
typedef signed int		S32;
typedef unsigned long long U64;
typedef unsigned long long ULL;
#endif



#endif /* TYPEDEFINITIONS_H_ */
