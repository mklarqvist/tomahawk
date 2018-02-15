#ifndef INDEX_INDEX_META_ENTRY_H_
#define INDEX_INDEX_META_ENTRY_H_

namespace Tomahawk{
namespace Totempole{

#define TWK_INDEX_META_ENTRY_SIZE (sizeof(U64)*4 + sizeof(U32)*2)

struct IndexMetaEntry{
public:
	typedef IndexMetaEntry self_type;

public:
	IndexMetaEntry() :
		index_begin(0),
		index_end(0),
		min_position(0),
		max_position(0),
		n_variants(0),
		uncompressed_size(0)
	{
	}

	IndexMetaEntry(const char* const data) :
		index_begin(*reinterpret_cast<const U64* const>(data)),
		index_end(*reinterpret_cast<const U64* const>(&data[sizeof(U32)])),
		min_position(*reinterpret_cast<const U64* const>(&data[sizeof(U32)*2+sizeof(U64)])),
		max_position(*reinterpret_cast<const U64* const>(&data[sizeof(U32)*2+sizeof(U64)])),
		n_variants(*reinterpret_cast<const U64* const>(&data[sizeof(U32)*2+sizeof(U64)*2])),
		uncompressed_size(*reinterpret_cast<const U64* const>(&data[sizeof(U32)*2+sizeof(U64)*3]))
	{
	}

	// Copy ctor
	IndexMetaEntry(const self_type& other) :
		index_begin(other.index_begin),
		index_end(other.index_end),
		min_position(other.min_position),
		max_position(other.max_position),
		n_variants(other.n_variants),
		uncompressed_size(other.uncompressed_size)
	{
	}
	~IndexMetaEntry() = default;

	inline U32 size(void) const{ return(this->n_variants); }
	inline const bool empty(void) const{ return(this->n_variants == 0); }
	inline void operator++(void){ ++this->n_variants; }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.index_begin << "-" << entry.index_end << '\t' << entry.min_position << '-' << entry.max_position << '\t' << entry.n_variants << '\t' << entry.uncompressed_size;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.index_begin),       sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.index_end),         sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.min_position),      sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.max_position),      sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),        sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U64));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.index_begin),       sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.index_end),         sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.min_position),      sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.max_position),      sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),        sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U64));
		return(stream);
	}

	void reset(void){
		this->index_begin       = 0;
		this->index_end         = 0;
		this->min_position      = 0;
		this->max_position      = 0;
		this->n_variants        = 0;
		this->uncompressed_size = 0;
	}

public:
	U32 index_begin;
	U32 index_end;
	U64 min_position;      // smallest bp position in tomahawk block
	U64 max_position;      // largest bp position in tomahawk block
	U64 n_variants;        // number of variants in this block
	U64 uncompressed_size; // uncompressed size of this block
};

}
}


#endif /* INDEX_INDEX_META_ENTRY_H_ */
