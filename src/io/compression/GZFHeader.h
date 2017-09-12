#ifndef TGZFHEADER_H_
#define TGZFHEADER_H_

#include "GZFConstants.h"

namespace Tomahawk{
namespace IO{

#pragma pack(1)
struct __headerBase{
private:
	typedef __headerBase self_type;

public:
	__headerBase();
	~__headerBase();

	BYTE ID1;
	BYTE ID2;
	BYTE CM;
	BYTE FLG;
	U32 MTIME;
	BYTE XFL;
	BYTE OS;
	U16 XLEN;
	BYTE SI1;
	BYTE SI2;
	U16 SLEN;

	friend std::ostream& operator<<(std::ostream& os, const self_type& header){
		os << "ID1\t" << (U32)header.ID1 << '\n';
		os << "ID2\t" << (U32)header.ID2 << '\n';
		os << "CM\t" << (U32)header.CM << '\n';
		os << "FLG\t" << (U32)header.FLG << '\n';
		os << "MTIME\t" << (S32)header.MTIME << '\n';
		os << "XFL\t" << (U32)header.XFL << '\n';
		os << "OS\t" << (U32)header.OS << '\n';
		os << "XLEN\t" << (U16)header.XLEN << '\n';
		os << "SI1\t" << (U32)header.SI1 << '\n';
		os << "SI2\t" << (U32)header.SI2 << '\n';
		os << "SLEN\t" << (U16)header.SLEN;

		return os;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		stream.read(reinterpret_cast<char*>(&header.ID1), sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&header.ID2), sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&header.CM), sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&header.FLG), sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&header.MTIME), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.XFL), sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&header.OS), sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&header.XLEN), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&header.SI1), sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&header.SI2), sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&header.SLEN), sizeof(U16));
		return(stream);
	}
};

/*
 TGZF header
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 84| 90|      4|        BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
  BGZF extension:
                ^                              ^   ^   ^
                |                              |   |   |
               FLG.EXTRA                     XLEN  T   Z
  TGZF format is compatible with GZIP. It limits the size of each compressed
  block to 2^32 bytes and adds and an extra "BC" field in the gzip header which
  records the size.
*/
#pragma pack(1)
struct TGZFHeader : public __headerBase{
private:
	typedef TGZFHeader self_type;
	typedef __headerBase parent_type;

public:
	U32 BSIZE;	// remainder size

	inline bool Validate(void) const{
		return(this->ID1 == Constants::GZIP_ID1
				&& this->ID2 == Constants::GZIP_ID2
				&& this->CM == Constants::CM_DEFLATE
				&& this->FLG == Constants::FLG_FEXTRA
				&& this->XLEN == Constants::TGZF_XLEN
				&& this->SI1 == Constants::TGZF_ID1
				&& this->SI2 == Constants::TGZF_ID2
				&& this->SLEN == Constants::TGZF_LEN
			);
	}

	friend std::ostream& operator<<(std::ostream& os, const self_type& header){
		const parent_type& b = *reinterpret_cast<const parent_type* const>(&header);
		os << b << '\n';
		os << "BSIZE\t" << (U32)header.BSIZE;
		return os;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		parent_type& b = *reinterpret_cast<parent_type*>(&header);
		stream >> b;
		stream.read(reinterpret_cast<char*>(&header.BSIZE), sizeof(U32));
		return(stream);
	}
};

/*
 BGZF/GZIP header (specialised from RFC 1952; little endian):
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 66| 67|      2|BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
  BGZF extension:
                ^                              ^   ^   ^
                |                              |   |   |
               FLG.EXTRA                     XLEN  B   C
  BGZF format is compatible with GZIP. It limits the size of each compressed
  block to 2^16 bytes and adds and an extra "BC" field in the gzip header which
  records the size.
*/
#pragma pack(1)
struct BGZFHeader : public __headerBase{
private:
	typedef BGZFHeader self_type;
	typedef __headerBase parent_type;

public:
	U16 BSIZE;	// remainder size

	inline bool Validate(void) const{
		return(this->ID1 == Constants::GZIP_ID1
				&& this->ID2 == Constants::GZIP_ID2
				&& this->CM == Constants::CM_DEFLATE
				&& this->FLG == Constants::FLG_FEXTRA
				&& this->XLEN == Constants::BGZF_XLEN
				&& this->SI1 == Constants::BGZF_ID1
				&& this->SI2 == Constants::BGZF_ID2
				&& this->SLEN == Constants::BGZF_LEN
			);
	}

	friend std::ostream& operator<<(std::ostream& os, const self_type& header){
		const parent_type& b = *reinterpret_cast<const parent_type* const>(&header);
		os << b << '\n';
		os << "BSIZE\t" << (U16)header.BSIZE;
		return os;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		parent_type& b = *reinterpret_cast<parent_type*>(&header);
		stream >> b;
		stream.read(reinterpret_cast<char*>(&header.BSIZE), sizeof(U16));
		return(stream);
	}
};

}
}


#endif /* TGZFHEADER_H_ */
