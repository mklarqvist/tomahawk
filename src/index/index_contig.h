#ifndef TOTEMPOLECONTIG_H_
#define TOTEMPOLECONTIG_H_

#include <ostream>
#include <fstream>

namespace Tomahawk{
namespace Totempole{

struct IndexContigBase{
public:
	typedef IndexContigBase self_type;

public:
	IndexContigBase(const U32& bases, const U32& n_char, const std::string& name) : n_bases(bases), n_char(n_char), name(name){}
	IndexContigBase() : n_bases(0), n_char(0){}
	~IndexContigBase(){}

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

public:
	U32 n_bases;      // length of contig
	U32 n_char;       // number of chars
	std::string name; // contig name
};

struct IndexContig : public IndexContigBase{
public:
	typedef IndexContig     self_type;
	typedef IndexContigBase parent_type;

public:
	IndexContig() : minPosition(0), maxPosition(0), blocksStart(0), blocksEnd(0){}
	~IndexContig(){}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.name << '\t' << entry.n_bases << '\t' << entry.minPosition << "-" << entry.maxPosition << '\t' << entry.blocksStart << "->" << entry.blocksEnd;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& base){
		const parent_type* const parent = reinterpret_cast<const parent_type* const>(&base);
		stream << *parent;

		stream.write(reinterpret_cast<const char*>(&base.minPosition), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.maxPosition), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.blocksStart), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.blocksEnd), sizeof(U32));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& base){
		parent_type* parent = reinterpret_cast<parent_type*>(&base);
		stream >> *parent;

		stream.read(reinterpret_cast<char*>(&base.minPosition), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&base.maxPosition), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&base.blocksStart), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&base.blocksEnd), sizeof(U32));
		return(stream);
	}

public:
	// contigID is implicit
	U32 minPosition;  // start position of contig
	U32 maxPosition;  // end position of contig
	U32 blocksStart;  // start IO-seek position of blocks
	U32 blocksEnd;    // end IO-seek position of blocks
};

}
}




#endif /* TOTEMPOLECONTIG_H_ */
