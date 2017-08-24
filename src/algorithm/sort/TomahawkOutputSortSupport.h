#ifndef TOMAHAWKOUTPUTSORTSUPPORT_H_
#define TOMAHAWKOUTPUTSORTSUPPORT_H_

namespace Tomahawk{
namespace IO {

struct PartialSortIndexHeader{
	typedef PartialSortIndexHeader self_type;

	friend std::istream& operator>>(std::istream& is, self_type& b){
		is.read(&b.header[0], Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH);
		char temp[sizeof(float)];
		is.read(&temp[0], sizeof(float));
		b.version = *reinterpret_cast<float*>(&temp[0]);
		return(is);
	}

	friend std::ofstream& operator<<(std::ofstream& of, self_type& b){
		of.write(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH);
		of.write((char*)&Tomahawk::Constants::PROGRAM_VERSION, sizeof(float));
		return(of);
	}

	bool validate(void) const{
		if(strncmp(this->header, Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH) != 0)
			return false;

		if(version < 0) return false;

		return true;
	}

	char header[Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH];
	float version;
};

#pragma pack(1)
struct PartialSortIndexHeaderEntry{
	typedef PartialSortIndexHeaderEntry self_type;

public:
	friend std::ostream& operator<<(std::ostream& os, const self_type& entry){
		os << entry.from << '\t' << entry.to;
		return(os);
	}

	friend std::ofstream& operator<<(std::ofstream& of, const self_type& entry){
		of.write((char*)&entry.from, sizeof(U64));
		of.write((char*)&entry.to, sizeof(U64));
		return(of);
	}

public:
	U64 from;
	U64 to;
};

}
}


#endif /* TOMAHAWKOUTPUTSORTSUPPORT_H_ */
