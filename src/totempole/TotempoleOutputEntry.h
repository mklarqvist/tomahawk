#ifndef TOTEMPOLEOUTPUTENTRY_H_
#define TOTEMPOLEOUTPUTENTRY_H_

#include "../support/type_definitions.h"

namespace Tomahawk {
namespace Totempole {

#pragma pack(push, 1)
struct __attribute__((packed, aligned(1))) TotempoleOutputEntryController{
	typedef TotempoleOutputEntryController self_type;

	BYTE sorted: 1,
	     expanded: 1,
         partial_sort: 1,
         unused: 5;

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write((const char*)reinterpret_cast<const char*>(&entry), 1);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry), 1);
		return(stream);
	}
};

struct __attribute__((packed, aligned(1))) TotempoleOutputEntry{
	typedef TotempoleOutputEntry self_type;

public:
	TotempoleOutputEntry() :
		byte_offset(0),
		byte_offset_end(0),
		n_entries(0),
		uncompressed_size(0)
	{}
	~TotempoleOutputEntry(){}

	inline const bool isValid(void) const{ return(this->byte_offset != 0); }
	inline void operator++(void){ ++this->n_entries; }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.byte_offset << '-' << entry.byte_offset_end << '\t' << entry.n_entries << '\t'
				<< entry.uncompressed_size;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.n_entries),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.n_entries),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));

		return(stream);
	}

	void reset(void){
		this->byte_offset       = 0;
		this->byte_offset_end   = 0;
		this->n_entries           = 0;
		this->uncompressed_size = 0;
	}

	// Do not return const_reference as it might be a temporary address
	inline const U32 size(void) const{ return(this->n_entries); }
	inline const U64 getStartOffset(void) const{ return(this->byte_offset); }
	inline const U64 getEndOffset(void) const{ return(this->byte_offset_end); }
	inline const U64 sizeCompressed(void) const{ return(this->byte_offset_end - this->byte_offset); }
	inline const U32 sizeUncompressed(void) const{ return(this->uncompressed_size); }
	inline void setStartOffset(const U64& byte_offset){ this->byte_offset = byte_offset; }
	inline void setEndOffset(const U64& byte_offset){ this->byte_offset_end = byte_offset_end; }
	inline void setSize(const U32& n_entries){ this->n_entries = n_entries; }
	inline void setUncompressedSize(const U32& l_bytes){ this->uncompressed_size = l_bytes; }

public:
	U64 byte_offset;		// tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;	// tellg() position in stream for start of record in Tomahawk file
	U32 n_entries; 			// number of variants in this block
	U32 uncompressed_size;	// uncompressed size of this block
};

struct __attribute__((packed, aligned(1))) TotempoleOutputSortedEntry{
	typedef TotempoleOutputSortedEntry self_type;

public:
	TotempoleOutputSortedEntry() :
		//contigID(0),
		fromBlock(-1),
		fromBlock_entries_offset(0),
		toBlock(-1),
		toBlock_entries_offset(0)
	{}
	~TotempoleOutputSortedEntry(){}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.fromBlock << ':' << entry.fromBlock_entries_offset << "->"
				<< entry.toBlock << ':' << entry.toBlock_entries_offset;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		//stream.write(reinterpret_cast<const char*>(&entry.contigID), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.fromBlock), sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.fromBlock_entries_offset), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.toBlock), sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.toBlock_entries_offset), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		//stream.read(reinterpret_cast<char*>(&entry.contigID), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.fromBlock), sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.fromBlock_entries_offset), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.toBlock), sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.toBlock_entries_offset), sizeof(U32));

		return(stream);
	}

	void reset(void){
		//this->contigID = 0;
		this->fromBlock                = -1;
		this->fromBlock_entries_offset = 0;
		this->toBlock                  = -1;
		this->toBlock_entries_offset   = 0;
	}

	inline void update(const U32& block, const U32& offset){
		this->toBlock = block;
		this->toBlock_entries_offset = offset;
	}

public:
	//U32 contigID;
	S32 fromBlock;	// tellg() position in stream for start of record in Tomahawk file
	U32 fromBlock_entries_offset;
	S32 toBlock; 	// tellg() position in stream for end of record in Tomahawk file
	U32 toBlock_entries_offset;
};
#pragma pack(pop)

}
}

#endif /* TOTEMPOLEOUTPUTENTRY_H_ */
