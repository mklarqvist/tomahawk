#ifndef MAGICCONSTANTS_H_
#define MAGICCONSTANTS_H_

#include <string>

#include "type_definitions.h"

extern int SILENT;

namespace Tomahawk{
namespace Constants{

extern std::string LITERAL_COMMAND_LINE;
extern std::string INTERPRETED_COMMAND;

// Versioning
const float PROGRAM_VERSION_MAJOR = 0.3; // major
const float PROGRAM_VERSION_MINOR = 3;

const double ALLOWED_ROUNDING_ERROR = 0.001;

const std::string PROGRAM_NAME = "tomahawk";
const std::string OUTPUT_SUFFIX = "twk";
const std::string OUTPUT_LD_SUFFIX = "two";

// Headers
const char* const WRITE_HEADER_MAGIC = "TOMAHAWK\1";
const U16 WRITE_HEADER_MAGIC_LENGTH = 9;

const BYTE TOMAHAWK_ALLELE_PACK_WIDTH = 2; // bit / allele
const BYTE TOMAHAWK_SNP_PACK_WIDTH = TOMAHAWK_ALLELE_PACK_WIDTH * 2; // bits / genotype

// Encoding for alleles
const char TOMAHAWK_ALLELE_LOOKUP[4] = {2, 3, 0, 1};
const char TOMAHAWK_ALLELE_LOOKUP_REVERSE[4] = {'0', '1', '.', '?'};

// Encoding for bases
const char* const REF_ALT_LOOKUP = "ATGCN";
const BYTE REF_ALT_A = 0;
const BYTE REF_ALT_T = 1;
const BYTE REF_ALT_G = 2;
const BYTE REF_ALT_C = 3;
const BYTE REF_ALT_N = 4;

// Upper bounds
// change to constants
const U32 UPPER_LIMIT_SAMPLES_8B =  ((1 << (8 - TOMAHAWK_SNP_PACK_WIDTH)) - 1); // 00001111 = 2^4 - 1
const U32 UPPER_LIMIT_SAMPLES_16B = ((1 << (16 - TOMAHAWK_SNP_PACK_WIDTH)) - 1); // 0000(1)12 = 2^12 - 1
const U32 UPPER_LIMIT_SAMPLES_32B = ((1 << (32 - TOMAHAWK_SNP_PACK_WIDTH)) - 1); // 0000(1)28 = 2^28 - 1
const U64 UPPER_LIMIT_SAMPLES_64B = (((U64)1 << (64 - TOMAHAWK_SNP_PACK_WIDTH)) - 1); // 0000(1)60 = 2^60 - 1

const BYTE SAMPLES_8B_MASK = 15;
const U16 SAMPLES_16B_MASK = 4095;
const U32 SAMPLES_32B_MASK = 268435455;
const U64 SAMPLES_64B_MASK = 1152921504606846976;

const BYTE eof_length = 64;
// EOF poem: "We will be known forever by the tracks we leave" - Santee Sioux Native Americans from Dakota;
const std::string eof_hex = "f3da5a14f8462d0e067eea643111437b3c033b61372ab0d55a45b5b1668f18db6a29d4c87b0c3ecdcaea374d936a406c248c851fe215c2c0669e2cfcd9f734a4";

}
}

#endif /* MAGICCONSTANTS_H_ */
