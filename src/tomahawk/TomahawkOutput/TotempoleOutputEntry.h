#ifndef TOTEMPOLEOUTPUTENTRY_H_
#define TOTEMPOLEOUTPUTENTRY_H_

namespace Tomahawk {
namespace Totempole {

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
		minPositionA(0),
		maxPositionA(0),
		contigIDB(-1),
		minPositionB(0),
		maxPositionB(0)
	{}
	~TotempoleOutputEntry(){}

	inline const bool isValid(void) const{ return(this->byte_offset != 0); }
	inline void operator++(void){ ++this->entries; }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.byte_offset << '-' << entry.byte_offset_end << '\t' << entry.entries << '\t'
				<< entry.uncompressed_size << '\t' << entry.contigIDA << '\t'
				<< entry.minPositionA << '-' << entry.maxPositionA << '\t'
				<< entry.contigIDB << '\t'
				<< entry.minPositionB << '-' << entry.maxPositionB;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.entries),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.contigIDA),    sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.minPositionA), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.maxPositionA), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.contigIDB),    sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.minPositionB), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.maxPositionB), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.entries),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.contigIDA),    sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.minPositionA), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.maxPositionA), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.contigIDB),    sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.minPositionB), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.maxPositionB), sizeof(U32));

		return(stream);
	}

	void reset(void){
		this->byte_offset = 0;
		this->byte_offset_end = 0;
		this->entries = 0;
		this->uncompressed_size = 0;
		this->contigIDA = -1;
		this->minPositionA = 0;
		this->maxPositionA = 0;
		this->contigIDB = -1;
		this->minPositionB = 0;
		this->maxPositionB = 0;
	}

public:
	U64 byte_offset;	// tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;	// tellg() position in stream for start of record in Tomahawk file
	U32 entries; 		// number of variants in this block
	U32 uncompressed_size; // uncompressed size of this block
	S32 contigIDA; // if contigID >= 0 then all entries belong to this contigID
	U32 minPositionA; // minPosition of entries. 0 if contigID = -1
	U32 maxPositionA; // maxPosition of entries, 0 if contigID = -1
	S32 contigIDB; // if contigID >= 0 then all entries belong to this contigID
	U32 minPositionB; // minPosition of entries. 0 if contigID = -1
	U32 maxPositionB; // maxPosition of entries, 0 if contigID = -1
};

}
}

#endif /* TOTEMPOLEOUTPUTENTRY_H_ */
