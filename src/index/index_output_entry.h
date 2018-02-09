#ifndef TOTEMPOLE_INDEX_OUTPUT_ENTRY_H_
#define TOTEMPOLE_INDEX_OUTPUT_ENTRY_H_

namespace Tomahawk{
namespace Totempole{

struct IndexOutputEntry{
	typedef IndexOutputEntry self_type;

public:
	IndexOutputEntry() :
		byte_offset(0),
		byte_offset_end(0),
		contigIDFrom(-1),
		contigIDTo(-1),
		minPositionFrom(-1),
		maxPositionFrom(-1),
		minPositionTo(-1),
		maxPositionTo(-1),
		n_entries(0),
		uncompressed_size(0)
	{

	}

	~IndexOutputEntry() = default;

	inline U32 size(void) const{ return(this->n_entries); }
	inline bool isValid(void) const{ return(this->byte_offset != 0); }
	inline void operator++(void){ ++this->n_entries; }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset),       sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end),   sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.contigIDFrom),      sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.contigIDTo),        sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.minPositionFrom),   sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.maxPositionFrom),   sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.minPositionTo),     sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.maxPositionTo),     sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.n_entries),         sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset),       sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end),   sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.contigIDFrom),      sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.contigIDTo),        sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.minPositionFrom),   sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.maxPositionFrom),   sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.minPositionTo),     sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.maxPositionTo),     sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.n_entries),         sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));
		return stream;
	}

	void reset(void){
		this->byte_offset       = 0;
		this->byte_offset_end   = 0;
		this->contigIDFrom      = 0;
		this->contigIDTo        = 0;
		this->minPositionFrom   = 0;
		this->maxPositionFrom   = 0;
		this->minPositionTo     = 0;
		this->maxPositionTo     = 0;
		this->n_entries         = 0;
		this->uncompressed_size = 0;
	}

public:
	U64 byte_offset;       // tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;   // tellg() position in stream for start of record in Tomahawk file
	S32 contigIDFrom;      // contig identifier
	S32 contigIDTo;        // contig identifier
	U64 minPositionFrom;   // smallest bp position in tomahawk block
	U64 maxPositionFrom;   // largest bp position in tomahawk block
	U64 minPositionTo;     // smallest bp position in tomahawk block
	U64 maxPositionTo;     // largest bp position in tomahawk block
	U32 n_entries;         // number of entries in this block
	U32 uncompressed_size; // uncompressed size of this block
};

}
}

#endif /* TOTEMPOLE_INDEX_OUTPUT_ENTRY_H_ */
