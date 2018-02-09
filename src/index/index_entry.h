#ifndef TOTEMPOLEENTRY_H_
#define TOTEMPOLEENTRY_H_

namespace Tomahawk{
namespace Totempole{

struct IndexEntry{
public:
	typedef IndexEntry self_type;

public:
	IndexEntry() : byte_offset(0), byte_offset_end(0), contigID(0), minPosition(0), maxPosition(0), n_variants(0), uncompressed_size(0){}

	// Copy ctor
	IndexEntry(const self_type& other) :
		byte_offset(other.byte_offset),
		byte_offset_end(other.byte_offset_end),
		contigID(other.contigID),
		minPosition(other.minPosition),
		maxPosition(other.maxPosition),
		n_variants(other.n_variants),
		uncompressed_size(other.uncompressed_size)
	{
	}
	~IndexEntry() = default;

	inline U32 size(void) const{ return(this->n_variants); }
	inline bool isValid(void) const{ return(this->byte_offset != 0); }
	inline void operator++(void){ ++this->n_variants; }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.byte_offset << '\t' << entry.byte_offset_end << '\t' << entry.contigID << '\t' << entry.minPosition << '-' << entry.maxPosition << '\t' << entry.n_variants << '\t' << entry.uncompressed_size;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset),       sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end),   sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.contigID),          sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition),       sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition),       sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),        sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset),       sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end),   sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.contigID),          sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.minPosition),       sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition),       sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),        sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));

		return(stream);
	}

	void reset(void){
		this->byte_offset       = 0;
		this->byte_offset_end   = 0;
		this->contigID          = 0;
		this->minPosition       = 0;
		this->maxPosition       = 0;
		this->n_variants        = 0;
		this->uncompressed_size = 0;
	}

public:
	U64 byte_offset;       // tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;   // tellg() position in stream for start of record in Tomahawk file
	S32 contigID;          // contig identifier
	U64 minPosition;       // smallest bp position in tomahawk block
	U64 maxPosition;       // largest bp position in tomahawk block
	U32 n_variants;        // number of variants in this block
	U32 uncompressed_size; // uncompressed size of this block
};

}
}

#endif /* TOTEMPOLEENTRY_H_ */
