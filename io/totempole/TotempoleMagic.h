#ifndef TOTEMPOLEMAGIC_H_
#define TOTEMPOLEMAGIC_H_

namespace Tomahawk{
namespace IO{

template <U16 length>
struct MAGICBase{
	typedef MAGICBase self_type;

public:
	MAGICBase(){}	// for reading
	MAGICBase(const char* target){ memcpy(&this->MAGIC[0], target, length); } // for writing
	virtual ~MAGICBase(){}

	friend std::istream& operator>>(std::istream& stream, self_type& base){
		stream.read(base.MAGIC, length);
		return(stream);
	}

	friend std::ostream& operator<<(std::ofstream& stream, const self_type& base){
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

	friend std::ostream& operator<<(std::ofstream& stream, const self_type& header){
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

	TomahawkOutputHeader() : n_contig(0){} // for reading
	TomahawkOutputHeader(const char* target, const U64 samples, const U32 n_contigs) :
		parent_type(target, samples),
		n_contig(n_contigs)
	{
		memcpy(&this->MAGIC[0], target, length);
	} // for writing

	friend std::ostream& operator<<(std::ofstream& stream, const self_type& header){
		stream.write(header.MAGIC, length);
		stream.write(reinterpret_cast<const char*>(&Tomahawk::Constants::PROGRAM_VERSION), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&header.samples), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&header.n_contig), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		stream.read(header.MAGIC, length);
		stream.read(reinterpret_cast<char *>(&header.version), sizeof(float));
		stream.read(reinterpret_cast<char *>(&header.samples), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&header.n_contig), sizeof(U32));
		return(stream);
	}

public:
	U32 n_contig;
};

}
}

#endif /* TOTEMPOLEMAGIC_H_ */
