#ifndef TOTEMPOLEHEADER_H_
#define TOTEMPOLEHEADER_H_

#include <bitset>

#include "support/magic_constants.h"

namespace tomahawk {
namespace totempole {

#define TWK_FOOTER_LENGTH (constants::eof_length + sizeof(U32) + sizeof(U64))

struct Footer{
public:
	typedef Footer self_type;

public:
	Footer() :
		offset_end_of_data(0),
		l_largest_uncompressed(0)
	{
		helpers::HexToBytes(constants::eof_hex, &this->EOF_marker[0]);
	}

	Footer(const char* const data) :
		offset_end_of_data(*reinterpret_cast<const U64* const>(data)),
		l_largest_uncompressed(*reinterpret_cast<const U64* const>(&data[sizeof(U64)]))
	{
		memcpy(&this->EOF_marker[0], &data[sizeof(U64)+sizeof(U32)], constants::eof_length);
	}

	~Footer() = default;

	inline const U64& getEODPosition(void) const{ return(this->offset_end_of_data); }
	inline const U32& getLargestUncompressedBlock(void) const{ return(this->l_largest_uncompressed); }
	inline U64& getEODPosition(void){ return(this->offset_end_of_data); }
	inline U32& getLargestUncompressedBlock(void){ return(this->l_largest_uncompressed); }

	inline const bool validate(void) const{
		if(this->offset_end_of_data == 0) return false;
		if(this->l_largest_uncompressed == 0) return false;

		BYTE reference[constants::eof_length];
		helpers::HexToBytes(constants::eof_hex, &reference[0]);

		if(strncmp(reinterpret_cast<const char* const>(&this->EOF_marker[0]), reinterpret_cast<const char* const>(&reference[0]), constants::eof_length) != 0) return false;
		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& footer){
		stream.write(reinterpret_cast<const char*>(&footer.offset_end_of_data),     sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&footer.l_largest_uncompressed), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&footer.EOF_marker[0]), constants::eof_length);
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& footer){
		stream.read(reinterpret_cast<char*>(&footer.offset_end_of_data),     sizeof(U64));
		stream.read(reinterpret_cast<char*>(&footer.l_largest_uncompressed), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&footer.EOF_marker[0]), constants::eof_length);
		return(stream);
	}

public:
	U64  offset_end_of_data;      // byte offset to end of Tomahawk data
	U32  l_largest_uncompressed;  // largest block-size in bytes
    BYTE EOF_marker[constants::eof_length];
};

}
}

#endif /* TOTEMPOLEHEADER_H_ */
