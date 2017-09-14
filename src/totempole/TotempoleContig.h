#ifndef TOTEMPOLECONTIG_H_
#define TOTEMPOLECONTIG_H_

#include <ostream>
#include <fstream>

namespace Tomahawk{
namespace Totempole{

struct TotempoleContigBase{
	typedef TotempoleContigBase self_type;

public:
	TotempoleContigBase(const U32& bases, const U32& n_char, const std::string& name) : bases(bases), n_char(n_char), name(name){}
	TotempoleContigBase() : bases(0), n_char(0){}
	~TotempoleContigBase(){}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.bases << '\t' << entry.n_char << '\t' << entry.name;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& base){
		stream.write(reinterpret_cast<const char*>(&base.bases), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.n_char), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&base.name[0]), base.name.size());
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& base){
		stream.read(reinterpret_cast<char *>(&base.bases), sizeof(U32));
		stream.read(reinterpret_cast<char *>(&base.n_char), sizeof(U32));
		base.name.resize(base.n_char);
		stream.read(&base.name[0], base.n_char);
		return(stream);
	}

public:
	U32 bases;			// length of contig
	U32 n_char;			// number of chars
	std::string name;	// contig name
};

struct TotempoleContig : public TotempoleContigBase{
typedef TotempoleContig self_type;

public:
	TotempoleContig() : minPosition(0), maxPosition(0), blocksStart(0), blocksEnd(0){}
	~TotempoleContig(){}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.name << '\t' << entry.bases << '\t' << entry.minPosition << "-" << entry.maxPosition << '\t' << entry.blocksStart << "->" << entry.blocksEnd;
		return stream;
	}

	// Updated second when read
	// contigID is implicit
	U32 minPosition;	// start position of contig
	U32 maxPosition;	// end position of contig
	U32 blocksStart; 	// start IO-seek position of blocks
	U32 blocksEnd;		// end IO-seek position of blocks
};

}
}




#endif /* TOTEMPOLECONTIG_H_ */
