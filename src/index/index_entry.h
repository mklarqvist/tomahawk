#ifndef TOTEMPOLEENTRY_H_
#define TOTEMPOLEENTRY_H_

#include <fstream>

#include "tomahawk/output_entry.h"

namespace tomahawk{
namespace totempole{

#define TWK_INDEX_ENTRY_SIZE (sizeof(U64)*4 + sizeof(S32) + sizeof(U32)*2)

struct IndexEntry{
public:
	typedef IndexEntry      self_type;
	typedef io::OutputEntry entry_type;

public:
	IndexEntry() :
		byte_offset(0),
		byte_offset_end(0),
		contigID(0),
		min_position(0),
		max_position(0),
		n_variants(0),
		uncompressed_size(0)
	{
	}

	IndexEntry(const char* const data) :
		byte_offset(*reinterpret_cast<const U64* const>(data)),
		byte_offset_end(*reinterpret_cast<const U64* const>(&data[sizeof(U64)])),
		contigID(*reinterpret_cast<const U64* const>(&data[sizeof(U64)*2])),
		min_position(*reinterpret_cast<const U64* const>(&data[sizeof(U64)*2+sizeof(S32)])),
		max_position(*reinterpret_cast<const U64* const>(&data[sizeof(U64)*2+sizeof(S32)+sizeof(U64)])),
		n_variants(*reinterpret_cast<const U64* const>(&data[sizeof(U64)*2+sizeof(S32)+sizeof(U64)*2])),
		uncompressed_size(*reinterpret_cast<const U64* const>(&data[sizeof(U64)*2+sizeof(S32)+sizeof(U64)*2+sizeof(U32)]))
	{
	}

	// Copy ctor
	IndexEntry(const self_type& other) :
		byte_offset(other.byte_offset),
		byte_offset_end(other.byte_offset_end),
		contigID(other.contigID),
		min_position(other.min_position),
		max_position(other.max_position),
		n_variants(other.n_variants),
		uncompressed_size(other.uncompressed_size)
	{
	}
	~IndexEntry() = default;

	inline U32 size(void) const{ return(this->n_variants); }
	inline bool isValid(void) const{ return(this->byte_offset != 0); }
	inline void operator++(void){ ++this->n_variants; }
	inline U64 sizeBytes(void) const{ return(this->byte_offset_end - this->byte_offset); }

	inline self_type& operator=(const entry_type& entry){
		this->contigID     = entry.AcontigID;
		this->min_position = entry.Aposition;
		this->max_position = entry.Aposition;
		return(*this);
	}

	void print(std::ostream& stream) const{
		stream << this->byte_offset << "->" << this->byte_offset_end <<
				" contig: " << this->contigID <<
				" min-max: " << this->min_position << "->" << this->max_position <<
				" variants: " << this->n_variants << " size: " << this->uncompressed_size;
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset),       sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end),   sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.contigID),          sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.min_position),      sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.max_position),      sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),        sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset),       sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end),   sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.contigID),          sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.min_position),      sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.max_position),      sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),        sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));

		return(stream);
	}

	void reset(void){
		this->byte_offset       = 0;
		this->byte_offset_end   = 0;
		this->contigID          = 0;
		this->min_position      = 0;
		this->max_position      = 0;
		this->n_variants        = 0;
		this->uncompressed_size = 0;
	}

	inline const bool overlaps(const S32& contigID) const{ return(this->contigID == contigID); }
	inline const bool overlaps(const S32& contigID, const U64& position) const{
		if(this->contigID != contigID) return false;
		return(position >= this->min_position && position <= this->max_position);
	}
	inline const bool overlaps(const S32& contigID, const U64& from_position, const U64& to_position) const{
		if(this->contigID != contigID) return false;
		if(to_position < this->min_position) return false;
		if(from_position > this->max_position) return false;
		return true;
	}

public:
	U64 byte_offset;       // tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;   // tellg() position in stream for start of record in Tomahawk file
	S32 contigID;          // contig identifier
	U64 min_position;      // smallest bp position in tomahawk block
	U64 max_position;      // largest bp position in tomahawk block
	U32 n_variants;        // number of variants in this block
	U32 uncompressed_size; // uncompressed size of this block
};

}
}

#endif /* TOTEMPOLEENTRY_H_ */
