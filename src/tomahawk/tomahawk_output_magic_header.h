#ifndef TOMAHAWK_TOMAHAWK_OUTPUT_MAGIC_HEADER_H_
#define TOMAHAWK_TOMAHAWK_OUTPUT_MAGIC_HEADER_H_

namespace Tomahawk{
namespace Base{

struct TomahawkOutputMagicHeader{
public:
	typedef TomahawkOutputMagicHeader self_type;

public:
	TomahawkOutputMagicHeader() :
		major_version(Tomahawk::Constants::PROGRAM_VERSION_MAJOR),
		minor_version(Tomahawk::Constants::PROGRAM_VERSION_MINOR),
		n_samples(0),
		n_contigs(0),
		controller(0)
	{
		memcpy(&this->magic_string[0],
               &Tomahawk::Constants::WRITE_HEADER_LD_MAGIC[0],
                Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH);
	}
	~TomahawkOutputMagicHeader() = default;

	inline bool validateMagic(void) const{ return(strncmp(&this->magic_string[0], &Tomahawk::Constants::WRITE_HEADER_LD_MAGIC[0], Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH) == 0); }
	inline bool validate(void) const{
		return(this->validateMagic() && this->n_samples > 0 && this->n_contigs > 0);
	}

	friend std::ostream& operator<<(std::ofstream& stream, const self_type& header){
		stream.write(header.magic_string, Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH);
		stream.write(reinterpret_cast<const char*>(&Tomahawk::Constants::PROGRAM_VERSION_MAJOR), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&Tomahawk::Constants::PROGRAM_VERSION_MINOR), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&header.n_samples), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&header.n_contigs), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.controller), sizeof(U16));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		stream.read(header.magic_string, Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH);
		stream.read(reinterpret_cast<char*>(&header.major_version), sizeof(float));
		stream.read(reinterpret_cast<char*>(&header.minor_version), sizeof(float));
		stream.read(reinterpret_cast<char*>(&header.n_samples), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&header.n_contigs), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.controller), sizeof(U16));
		return(stream);
	}

public:
	char magic_string[Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH];
	float major_version;
	float minor_version;
	U64 n_samples;
	U32 n_contigs;
	U16 controller;
};

}
}



#endif /* TOMAHAWK_TOMAHAWK_OUTPUT_MAGIC_HEADER_H_ */
