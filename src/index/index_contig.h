#ifndef TOTEMPOLECONTIG_H_
#define TOTEMPOLECONTIG_H_

#include <ostream>
#include <fstream>

#include "io/basic_buffer.h"
#include "support/type_definitions.h"

namespace tomahawk{
namespace totempole{

struct HeaderContig{
public:
	typedef HeaderContig    self_type;
	typedef io::BasicBuffer buffer_type;

public:
	HeaderContig(const U32& bases, const U32& n_char, const std::string& name) :
		n_bases(bases),
		n_char(n_char),
		name(name)
	{}

	HeaderContig() : n_bases(0), n_char(0){}

	HeaderContig(const char* const data) :
		n_bases(*reinterpret_cast<const U32* const>(data)),
		n_char(*reinterpret_cast<const U32* const>(&data[sizeof(U32)]))
	{
		this->name.resize(this->n_char);
		memcpy(&this->name[0], &data[sizeof(U32)+sizeof(U32)], this->n_char);
	}

	~HeaderContig(){}

	const U32 interpret(const char* const data){
		this->n_bases = *reinterpret_cast<const U32* const>(data);
		this->n_char = *reinterpret_cast<const U32* const>(&data[sizeof(U32)]);
		this->name.resize(this->n_char);
		memcpy(&this->name[0], &data[sizeof(U32)+sizeof(U32)], this->n_char);
		return(sizeof(U32) + sizeof(U32) + this->n_char);
	}

	const U32 interpret(const U32& bases, const U32& n_char, const std::string& name){
		this->n_bases = bases;
		this->n_char  = n_char;
		this->name    = name;
		return(sizeof(U32) + sizeof(U32) + this->n_char);
	}

	inline const bool operator==(const self_type& other) const{
		if(this->n_bases != other.n_bases) return false;
		if(this->n_char != other.n_char) return false;
		if(this->name != other.name) return false;
		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.n_bases << '\t' << entry.n_char << '\t' << entry.name;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& base){
		stream.write(reinterpret_cast<const char*>(&base.n_bases), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.n_char), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.name[0]), base.name.size());
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& base){
		stream.read(reinterpret_cast<char *>(&base.n_bases), sizeof(U32));
		stream.read(reinterpret_cast<char *>(&base.n_char), sizeof(U32));
		base.name.resize(base.n_char);
		stream.read(&base.name[0], base.n_char);
		return(stream);
	}

	friend buffer_type& operator+=(buffer_type& buffer, self_type& base){
		buffer += base.n_bases;
		buffer += base.n_char;
		buffer.Add(base.name.data(), base.name.size());
		return(buffer);
	}

public:
	U32 n_bases;      // length of contig
	U32 n_char;       // number of chars
	std::string name; // contig name
};

struct IndexContig : public HeaderContig{
public:
	typedef IndexContig     self_type;
	typedef HeaderContig    parent_type;

public:
	IndexContig() : min_position(0), max_position(0), blocks_start(0), blocks_end(0){}
	~IndexContig(){}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.name << '\t' << entry.n_bases << '\t' << entry.min_position << "-" << entry.max_position << '\t' << entry.blocks_start << "->" << entry.blocks_end;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& base){
		const parent_type* const parent = reinterpret_cast<const parent_type* const>(&base);
		stream << *parent;

		stream.write(reinterpret_cast<const char*>(&base.min_position), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.max_position), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.blocks_start), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.blocks_end), sizeof(U32));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& base){
		parent_type* parent = reinterpret_cast<parent_type*>(&base);
		stream >> *parent;

		stream.read(reinterpret_cast<char*>(&base.min_position), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&base.max_position), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&base.blocks_start), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&base.blocks_end), sizeof(U32));
		return(stream);
	}

public:
	// contigID is implicit
	U32 min_position;  // start position of contig
	U32 max_position;  // end position of contig
	U32 blocks_start;  // start IO-seek position of blocks
	U32 blocks_end;    // end IO-seek position of blocks
};

}
}




#endif /* TOTEMPOLECONTIG_H_ */
