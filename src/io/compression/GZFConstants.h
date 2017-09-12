#ifndef GZFCONSTANTS_H_
#define GZFCONSTANTS_H_

#include "../../support/TypeDefinitions.h"

namespace Tomahawk{
namespace IO{
namespace Constants{

// zlib & TGZF constants
const BYTE GZIP_ID1   = 31;
const BYTE GZIP_ID2   = 139;
const BYTE CM_DEFLATE = 8;
const BYTE FLG_FEXTRA = 4;
const BYTE OS_UNKNOWN = 255;
const BYTE TGZF_XLEN  = 8;
const BYTE TGZF_ID1   = 84;
const BYTE TGZF_ID2   = 90;
const BYTE TGZF_LEN   = 4;
const BYTE BGZF_XLEN  = 6;
const BYTE BGZF_ID1   = 66;
const BYTE BGZF_ID2   = 67;
const BYTE BGZF_LEN   = 2;

const SBYTE		GZIP_WINDOW_BITS			= -15;
const SBYTE		Z_DEFAULT_MEM_LEVEL			= 8;
const BYTE  	TGZF_BLOCK_HEADER_LENGTH  	= 20;
const BYTE		TGZF_BLOCK_FOOTER_LENGTH  	= 8;
const BYTE  	BGZF_BLOCK_HEADER_LENGTH  	= 18;

}
}
}

#endif /* GZFCONSTANTS_H_ */
