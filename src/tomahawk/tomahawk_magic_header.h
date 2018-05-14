#ifndef TOMAHAWK_TOMAHAWK_MAGIC_HEADER_H_
#define TOMAHAWK_TOMAHAWK_MAGIC_HEADER_H_

#include "../support/MagicConstants.h"

namespace tomahawk{
namespace base{

struct TomahawkMagicHeader{
public:
	typedef TomahawkMagicHeader self_type;

public:
	TomahawkMagicHeader() :
		major_version(constants::PROGRAM_VERSION_MAJOR),
		minor_version(constants::PROGRAM_VERSION_MINOR),
		file_type(0),
		n_samples(0),
		n_contigs(0),
		controller(0),
		l_header(0),
		l_header_uncompressed(0)
	{
		memcpy(&this->magic_string[0],
               &constants::WRITE_HEADER_MAGIC[0],
                constants::WRITE_HEADER_MAGIC_LENGTH);
	}

	TomahawkMagicHeader(const self_type& other) :
		major_version(other.major_version),
		minor_version(other.minor_version),
		file_type(other.file_type),
		n_samples(other.n_samples),
		n_contigs(other.n_contigs),
		controller(other.controller),
		l_header(0),
		l_header_uncompressed(0)
	{
		memcpy(&this->magic_string[0],
			   &other.magic_string[0],
				constants::WRITE_HEADER_MAGIC_LENGTH);
	}

	~TomahawkMagicHeader() = default;

	inline const U64& getNumberSamples(void) const{ return(this->n_samples); }
	inline U64& getNumberSamples(void){ return(this->n_samples); }
	inline const U32& getNumberContigs(void) const{ return(this->n_contigs); }
	inline U32& getNumberContigs(void){ return(this->n_contigs); }

	inline bool validateMagic(void) const{ return(strncmp(&this->magic_string[0], &constants::WRITE_HEADER_MAGIC[0], constants::WRITE_HEADER_MAGIC_LENGTH) == 0); }
	inline bool validate(void) const{
		return(this->validateMagic() && this->n_samples > 0 && this->n_contigs > 0 && (this->major_version > 0 || this->minor_version > 0) && this->l_header > 0 && this->l_header_uncompressed > 0);
	}

	inline const bool operator==(const self_type& other) const{
		if(strncmp(&this->magic_string[0], &other.magic_string[0], constants::WRITE_HEADER_MAGIC_LENGTH) != 0) return false;
		if(this->file_type != other.file_type) return false;
		if(this->n_samples != other.n_samples) return false;
		if(this->n_contigs != other.n_contigs) return false;
		return true;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& header){
		stream.write(header.magic_string, constants::WRITE_HEADER_MAGIC_LENGTH);
		stream.write(reinterpret_cast<const char*>(&constants::PROGRAM_VERSION_MAJOR), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&constants::PROGRAM_VERSION_MINOR), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&header.file_type),  sizeof(BYTE));
		stream.write(reinterpret_cast<const char*>(&header.n_samples),  sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&header.n_contigs),  sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.controller), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&header.l_header),   sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.l_header_uncompressed), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		stream.read(header.magic_string, constants::WRITE_HEADER_MAGIC_LENGTH);
		stream.read(reinterpret_cast<char*>(&header.major_version), sizeof(float));
		stream.read(reinterpret_cast<char*>(&header.minor_version), sizeof(float));
		stream.read(reinterpret_cast<char*>(&header.file_type),     sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&header.n_samples),     sizeof(U64));
		stream.read(reinterpret_cast<char*>(&header.n_contigs),     sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.controller),    sizeof(U16));
		stream.read(reinterpret_cast<char*>(&header.l_header),      sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.l_header_uncompressed), sizeof(U32));
		return(stream);
	}

public:
	char  magic_string[constants::WRITE_HEADER_MAGIC_LENGTH];
	float major_version;
	float minor_version;
	BYTE  file_type;
	U64   n_samples;
	U32   n_contigs;
	U16   controller;
	U32   l_header;
	U32   l_header_uncompressed;
};

}
}

#endif /* TOMAHAWK_TOMAHAWK_MAGIC_HEADER_H_ */
