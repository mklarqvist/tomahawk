#ifndef TOTEMPOLEMAGIC_H_
#define TOTEMPOLEMAGIC_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include "TotempoleOutputEntry.h"
#include "../support/MagicConstants.h"

namespace Tomahawk{
namespace IO{

template <U16 length>
struct MAGICBase{
	typedef MAGICBase self_type;

public:
	MAGICBase(){}	// for reading
	MAGICBase(const char* target){ memcpy(&this->MAGIC[0], target, length); } // for writing
	MAGICBase(const self_type& other){ memcpy(&this->MAGIC[0], &other.MAGIC[0], length); }
	virtual ~MAGICBase(){}

	friend std::istream& operator>>(std::istream& stream, self_type& base){
		stream.read(base.MAGIC, length);
		return(stream);
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& base){
		stream.write(base.MAGIC, length);
		return stream;
	}

	virtual inline bool validate(const char* match) const{ return(strncmp(&this->MAGIC[0], match, length) == 0); }

public:
	char MAGIC[length];
};

template <U16 length>
struct TomahawkHeader : public MAGICBase<length>{
	typedef TomahawkHeader self_type;
	typedef MAGICBase<length> parent_type;

	TomahawkHeader() : version(0), samples(0){} // for reading
	TomahawkHeader(const char* target, const U64 samples) :
		version(Tomahawk::Constants::PROGRAM_VERSION),
		samples(samples)
	{
		memcpy(&this->MAGIC[0], target, length);
	} // for writing

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& header){
		stream.write(header.MAGIC, length);
		stream.write(reinterpret_cast<const char*>(&Tomahawk::Constants::PROGRAM_VERSION), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&header.samples), sizeof(U64));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		stream.read(header.MAGIC, length);
		stream.read(reinterpret_cast<char *>(&header.version), sizeof(float));
		stream.read(reinterpret_cast<char *>(&header.samples), sizeof(U64));
		return(stream);
	}

public:
	float version;
	U64 samples;
};

template <U16 length>
struct TomahawkOutputHeader : public TomahawkHeader<length>{
	typedef TomahawkOutputHeader self_type;
	typedef TomahawkHeader<length> parent_type;
	typedef Totempole::TotempoleOutputEntryController totempole_controller_byte;

	TomahawkOutputHeader() : n_contig(0), n_entries(0){} // for reading
	TomahawkOutputHeader(const char* target, const U64 samples, const U32 n_contigs) :
		parent_type(target, samples),
		n_contig(n_contigs),
		n_entries(0)
	{
		memcpy(&this->MAGIC[0], target, length);
	} // for writing

	TomahawkOutputHeader(const self_type& other) :
		parent_type(other.MAGIC, other.samples),
		n_contig(other.n_contig),
		n_entries(other.n_entries),
		controller(other.controller)
	{

	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& header){
		stream.write(header.MAGIC, length);
		stream.write(reinterpret_cast<const char*>(&Tomahawk::Constants::PROGRAM_VERSION), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&header.samples), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&header.n_contig), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.n_entries), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.controller), sizeof(BYTE));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		stream.read(header.MAGIC, length);
		stream.read(reinterpret_cast<char *>(&header.version), sizeof(float));
		stream.read(reinterpret_cast<char *>(&header.samples), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&header.n_contig), sizeof(U32));

		// Legacy fix
		if(header.version >= 0.2){
			stream.read(reinterpret_cast<char*>(&header.n_entries), sizeof(U32));
			stream.read(reinterpret_cast<char*>(&header.controller), sizeof(BYTE));
		} else
			header.n_entries = 0;

		return(stream);
	}

public:
	U32 n_contig;
	U32 n_entries;
	totempole_controller_byte controller;
};

}
}

#endif /* TOTEMPOLEMAGIC_H_ */
