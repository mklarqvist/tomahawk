#ifndef TOTEMPOLEOUTPUTENTRY_H_
#define TOTEMPOLEOUTPUTENTRY_H_

#include "../support/TypeDefinitions.h"

namespace Tomahawk {
namespace Totempole {

#pragma pack(1)
struct TotempoleOutputEntryController{
private:
	typedef TotempoleOutputEntryController self_type;

public:
	TotempoleOutputEntryController(void) : sorted(0), expanded(0), partial_sort(0), unused(0){}
	~TotempoleOutputEntryController(){}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write((const char*)reinterpret_cast<const char*>(&entry), 1);
		return(stream);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << "sorted=" << (int)entry.sorted << "; expanded=" << (int)entry.expanded << ";partial=" << (int)entry.partial_sort;
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry), 1);
		return(stream);
	}

public:
	BYTE sorted: 1,
	     expanded: 1,
		 partial_sort: 1,
		 unused: 5;
};

#pragma pack(1)
struct TotempoleOutputEntry{
	typedef TotempoleOutputEntry self_type;

public:
	TotempoleOutputEntry() :
		byte_offset(0),
		byte_offset_end(0),
		entries(0),
		uncompressed_size(0)
	{}
	~TotempoleOutputEntry(){}

	inline const bool isValid(void) const{ return(this->byte_offset != 0); }
	inline void operator++(void){ ++this->entries; }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.byte_offset << '-' << entry.byte_offset_end << '\t' << entry.entries << '\t'
				<< entry.uncompressed_size;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.entries),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.entries),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));

		return(stream);
	}

	void reset(void){
		this->byte_offset = 0;
		this->byte_offset_end = 0;
		this->entries = 0;
		this->uncompressed_size = 0;
	}

public:
	U64 byte_offset;		// tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;	// tellg() position in stream for start of record in Tomahawk file
	U32 entries; 			// number of variants in this block
	U32 uncompressed_size;	// uncompressed size of this block
};

#pragma pack(1)
struct TotempoleOutputSortedEntry{
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
		this->fromBlock = -1;
		this->fromBlock_entries_offset = 0;
		this->toBlock = -1;
		this->toBlock_entries_offset = 0;
	}

	inline void update(const U32& block, const U32& offset){
		this->toBlock = block;
		this->toBlock_entries_offset = offset;
	}

public:
	//U32 contigID;	// tellg() position in stream for start of record in Tomahawk file
	S32 fromBlock;	// tellg() position in stream for start of record in Tomahawk file
	U32 fromBlock_entries_offset;
	S32 toBlock; 	// number of variants in this block
	U32 toBlock_entries_offset;
};

}
}

#endif /* TOTEMPOLEOUTPUTENTRY_H_ */
