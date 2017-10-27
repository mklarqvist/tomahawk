#ifndef MAGICCONSTANTS_H_
#define MAGICCONSTANTS_H_

#include <string>
#include "../support/TypeDefinitions.h"

extern int SILENT;

namespace Tomahawk{
namespace Constants{

extern std::string LITERAL_COMMAND_LINE;
extern std::string INTERPRETED_COMMAND;

// Versioning
const float PROGRAM_VERSION = 0.4; // major

const double ALLOWED_ROUNDING_ERROR = 0.001;

const std::string PROGRAM_NAME = "tomahawk";
const std::string OUTPUT_SUFFIX = "twk";
const std::string OUTPUT_INDEX_SUFFIX = "twi";
const std::string OUTPUT_LD_SUFFIX = "two";
const std::string OUTPUT_LD_PARTIAL_SORT_INDEX_SUFFIX = "twsi";
const std::string OUTPUT_LD_SORT_INDEX_SUFFIX = "toi";

// Headers
const char* const WRITE_HEADER_MAGIC = "TOMAHAWK\1";
const char* const WRITE_HEADER_INDEX_MAGIC = "TOTEMPOLE\1";
const char* const WRITE_HEADER_LD_MAGIC = "TOMAHAWK~OUTPUT\1";
const char* const WRITE_HEADER_LD_SORT_MAGIC = "TOMAHAWK~OUTPUT~INDEX\1";
const U16 WRITE_HEADER_MAGIC_LENGTH = 9;
const U16 WRITE_HEADER_MAGIC_INDEX_LENGTH = 10;
const U16 WRITE_HEADER_LD_MAGIC_LENGTH = 16;
const U16 WRITE_HEADER_LD_SORT_MAGIC_LENGTH = 22;

const BYTE TOMAHAWK_ALLELE_PACK_WIDTH = 2; // bit / allele
const BYTE TOMAHAWK_SNP_PACK_WIDTH = TOMAHAWK_ALLELE_PACK_WIDTH * 2; // bits / genotype
const BYTE TOMAHAWK_SHIFT_SIZE = Constants::TOMAHAWK_SNP_PACK_WIDTH + 1;

// Encoding for alleles
const char TOMAHAWK_ALLELE_LOOKUP[4] = {2, 3, 0, 1};
const char TOMAHAWK_ALLELE_LOOKUP_REVERSE[4] = {'0', '1', '.', '?'};

// 0 -> 0
// 1 -> 1
// 2 -> 4
// 3 -> 5
const BYTE TOMAHAWK_ALLELE_REDUCED_MAP[4] = {0, 1, 4, 5};
// 0 -> 0
// 1 -> 1
// 4 -> 4
// 5 -> 5
const BYTE TOMAHAWK_ALLELE_SELF_MAP[6] = {0, 1, 0, 0, 4, 5};

// Encoding for bases
const char* const REF_ALT_LOOKUP = "ATGCN";
const BYTE REF_ALT_A = 0;
const BYTE REF_ALT_T = 1;
const BYTE REF_ALT_G = 2;
const BYTE REF_ALT_C = 3;
const BYTE REF_ALT_N = 4;

// Upper bounds
// change to constants
const U32 UPPER_LIMIT_SAMPLES_8B = ((1 << (8 - (TOMAHAWK_SNP_PACK_WIDTH + 1))) - 1);   // 00001111 = 2^4 - 1
const U32 UPPER_LIMIT_SAMPLES_16B = ((1 << (16 - (TOMAHAWK_SNP_PACK_WIDTH + 1))) - 1); // 0000(1)12 = 2^12 - 1
const U32 UPPER_LIMIT_SAMPLES_32B = ((1 << (32 - (TOMAHAWK_SNP_PACK_WIDTH + 1))) - 1); // 0000(1)28 = 2^28 - 1
const U64 UPPER_LIMIT_SAMPLES_64B = (((U64)1 << (64 - (TOMAHAWK_SNP_PACK_WIDTH + 1))) - 1); // 0000(1)60 = 2^60 - 1

const BYTE SAMPLES_8B_MASK = 7;
const U16 SAMPLES_16B_MASK = 2047;
const U32 SAMPLES_32B_MASK = 134217727;
const U64 SAMPLES_64B_MASK = 576460752303423488;

// EOF
//const char* const TOMAHAWK_EOF_MARKER = "We will be known forever by the tracks we leave" - Santee Sioux Native Americans from Dakota;
//const U32 TOMAHAWK_EOF_MARKER_LENGTH = 31;

const BYTE eof_length = 6;
const U64 eof[6] = {2336361506924422487, 7959953386435011938, 8243124871055238688, 2334386829831791136, 8583987794834190964, 28464622577219173};
// EOF poem: "We will be known forever by the tracks we leave"

}
}

#endif /* MAGICCONSTANTS_H_ */
