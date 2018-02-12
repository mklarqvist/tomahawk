#ifndef TOTEMPOLE_INDEX_OUTPUT_ENTRY_H_
#define TOTEMPOLE_INDEX_OUTPUT_ENTRY_H_

namespace Tomahawk{
namespace Totempole{

struct IndexOutputEntry{
	typedef IndexOutputEntry self_type;

public:
	IndexOutputEntry() :
		byte_offset_from(0),
		byte_offset_to(0),
		contigID_from(-1),
		contigID_to(-1),
		min_position_from(-1),
		max_position_from(-1),
		min_position_to(-1),
		max_position_to(-1),
		n_entries(0),
		uncompressed_size(0)
	{

	}

	~IndexOutputEntry() = default;

	inline const U32& size(void) const{ return(this->n_entries); }
	inline bool isValid(void) const{ return(this->byte_offset_from != 0); }
	inline void operator++(void){ ++this->n_entries; }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_from),  sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_to),    sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.contigID_from),     sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.contigID_to),       sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.min_position_from), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.max_position_from), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.min_position_to),   sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.max_position_to),   sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.n_entries),         sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_from),  sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_to),    sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.contigID_from),     sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.contigID_to),       sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.min_position_from), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.max_position_from), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.min_position_to),   sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.max_position_to),   sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.n_entries),         sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));
		return stream;
	}

	void reset(void){
		this->byte_offset_from       = 0;
		this->byte_offset_to   = 0;
		this->contigID_from      = 0;
		this->contigID_to        = 0;
		this->min_position_from   = 0;
		this->max_position_from   = 0;
		this->min_position_to     = 0;
		this->max_position_to     = 0;
		this->n_entries         = 0;
		this->uncompressed_size = 0;
	}

public:
	U64 byte_offset_from;       // tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_to;   // tellg() position in stream for start of record in Tomahawk file
	S32 contigID_from;      // contig identifier
	S32 contigID_to;        // contig identifier
	U64 min_position_from;   // smallest bp position in tomahawk block
	U64 max_position_from;   // largest bp position in tomahawk block
	U64 min_position_to;     // smallest bp position in tomahawk block
	U64 max_position_to;     // largest bp position in tomahawk block
	U32 n_entries;         // number of entries in this block
	U32 uncompressed_size; // uncompressed size of this block
};

}
}

#endif /* TOTEMPOLE_INDEX_OUTPUT_ENTRY_H_ */
