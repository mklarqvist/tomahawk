#ifndef TOTEMPOLEENTRY_H_
#define TOTEMPOLEENTRY_H_

namespace Tomahawk{
namespace Totempole{

#pragma pack(push, 1)
struct __attribute__((packed, aligned(1))) TotempoleEntry{
	typedef TotempoleEntry self_type;

public:
	TotempoleEntry() : byte_offset(0), byte_offset_end(0), contigID(0), minPosition(0), maxPosition(0), variants(0), uncompressed_size(0){}
	TotempoleEntry(const self_type& other) :
		byte_offset(other.byte_offset),
		byte_offset_end(other.byte_offset_end),
		contigID(other.contigID),
		minPosition(other.minPosition),
		maxPosition(other.maxPosition),
		variants(other.variants),
		uncompressed_size(other.uncompressed_size)
	{
	}
	~TotempoleEntry(){}

	inline bool isValid(void) const{ return(this->byte_offset != 0); }
	inline void operator++(void){ ++this->variants; }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.byte_offset << '\t' << entry.byte_offset_end << '\t' << entry.contigID << '\t' << entry.minPosition << '-' << entry.maxPosition << '\t' << entry.variants << '\t' << entry.uncompressed_size;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.contigID),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.variants),    sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.contigID),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.minPosition), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.variants),    sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));

		return(stream);
	}

	void reset(void){
		this->byte_offset = 0;
		this->byte_offset_end = 0;
		this->contigID = 0;
		this->minPosition = 0;
		this->maxPosition = 0;
		this->variants = 0;
		this->uncompressed_size = 0;
	}

public:
	U64 byte_offset;	// tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;// tellg() position in stream for start of record in Tomahawk file
	S32 contigID;		// contig identifier
	U32 minPosition;	// smallest bp position in tomahawk block
	U32 maxPosition;	// largest bp position in tomahawk block
	U16 variants; 		// number of variants in this block
	U32 uncompressed_size; // uncompressed size of this block
};
#pragma pack(pop)

}
}

#endif /* TOTEMPOLEENTRY_H_ */
