#ifndef TOTEMPOLEOUTPUTENTRY_H_
#define TOTEMPOLEOUTPUTENTRY_H_

#include "../support/TypeDefinitions.h"

namespace Tomahawk {
namespace Totempole {

#pragma pack(1)
struct TotempoleOutputEntryController{
	typedef TotempoleOutputEntryController self_type;

	BYTE sorted: 1, unused: 7;

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write((const char*)reinterpret_cast<const char*>(&entry), 1);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry), 1);
		return(stream);
	}
};

#pragma pack(1)
struct TotempoleOutputEntry{
	typedef TotempoleOutputEntry self_type;

public:
	TotempoleOutputEntry() :
		byte_offset(0),
		byte_offset_end(0),
		entries(0),
		uncompressed_size(0),
		contigIDA(-1),
		minPositionA(-1),
		maxPositionA(-1),
		contigIDB(-1),
		minPositionB(-1),
		maxPositionB(-1)
	{}
	~TotempoleOutputEntry(){}

	inline const bool isValid(void) const{ return(this->byte_offset != 0); }
	inline void operator++(void){ ++this->entries; }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.byte_offset << '-' << entry.byte_offset_end << '\t' << entry.entries << '\t'
				<< entry.uncompressed_size << '\t'
				<< entry.contigIDA << '\t' << entry.minPositionA << '-' << entry.maxPositionA << '\t'
				<< entry.contigIDB << '\t' << entry.minPositionB << '-' << entry.maxPositionB;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.entries),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.contigIDA),    sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.minPositionA), sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.maxPositionA), sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.contigIDB),    sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.minPositionB), sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.maxPositionB), sizeof(S32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.entries),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.contigIDA),    sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.minPositionA), sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.maxPositionA), sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.contigIDB),    sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.minPositionB), sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.maxPositionB), sizeof(S32));

		return(stream);
	}

	void reset(void){
		this->byte_offset = 0;
		this->byte_offset_end = 0;
		this->entries = 0;
		this->uncompressed_size = 0;
		this->contigIDA = -1;
		this->minPositionA = -1;
		this->maxPositionA = -1;
		this->contigIDB = -1;
		this->minPositionB = -1;
		this->maxPositionB = -1;
	}

public:
	U64 byte_offset;		// tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;	// tellg() position in stream for start of record in Tomahawk file
	U32 entries; 			// number of variants in this block
	U32 uncompressed_size;	// uncompressed size of this block
	S32 contigIDA; 			// if contigID >= 0 then all entries belong to this contigID
	S32 minPositionA; 		// minPosition of entries. 0 if contigID = -1 or not sorted
	S32 maxPositionA;		// maxPosition of entries, 0 if contigID = -1 or not sorted
	S32 contigIDB;			// if contigID >= 0 then all entries belong to this contigID
	S32 minPositionB;		// minPosition of entries. 0 if contigID = -1 or not sorted
	S32 maxPositionB;		// maxPosition of entries, 0 if contigID = -1 or not sorted
};

}
}

#endif /* TOTEMPOLEOUTPUTENTRY_H_ */
