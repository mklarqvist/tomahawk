#ifndef VCFHEADERCONSTANTS_H_
#define VCFHEADERCONSTANTS_H_

#include <string>
#include "../../support/TypeDefinitions.h"

namespace Tomahawk{
namespace VCF{
namespace Constants{

const std::string HEADER_COLUMN = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
const char VCF_DELIMITER = '\t';

const std::string HEADER_VCF_FORMAT = "##fileformat=VCFv";
const std::string HEADER_VCF_VERSION = "##fileformat=VCFv4";

enum HEADER_BASE_TYPE { Fileformat, Contig, Alt, Info, Filter, Format };
const std::string HEADER_FILEFORMAT = "##fileformat=";
const std::string HEADER_CONTIG = "##contig=";
const std::string HEADER_ALT = "##alt=";
const std::string HEADER_INFO = "##info=";
const std::string HEADER_FILTER = "##filter=";
const std::string HEADER_FORMAT = "##format=";

enum HEADER_TYPE { Integer, Float, Character, String, Flag };
const std::string HEADER_TYPE_INTEGER = ",Type=Integer,";
const std::string HEADER_TYPE_FLOAT = ",Type=Float,";
const std::string HEADER_TYPE_CHARACTER = ",Type=Character,";
const std::string HEADER_TYPE_STRING = ",Type=String,";
const std::string HEADER_TYPE_FLAG = ",Type=Flag,";
const std::string HEADER_FIELD_NUMBER = "Number=";
const std::string HEADER_FIELD_DESCRIPTION = "Description=\"";

// Case sensitive
const std::string HEADER_ALT_TYPE_DEL = "ID=DEL,";
const std::string HEADER_ALT_TYPE_INS = "ID=INS,";
const std::string HEADER_ALT_TYPE_DUP = "ID=DUP,";
const std::string HEADER_ALT_TYPE_INV = "ID=INV,";
const std::string HEADER_ALT_TYPE_CNV = "ID=CNV,";
const std::string HEADER_ALT_TYPE_DUP_TANDEM = "ID=DUP:TANDEM,";
const std::string HEADER_ALT_TYPE_DUP_ME = "ID=DUP:ME,";
const std::string HEADER_ALT_TYPE_INS_ME = "ID=INS:ME,";

// GT
const std::string GT_ONLY = "GT";

}
}
}


#endif /* VCFHEADERCONSTANTS_H_ */
