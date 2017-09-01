#ifndef TOTEMPOLEOUTPUTENTRY_H_
#define TOTEMPOLEOUTPUTENTRY_H_

#pragma pack(1)
struct TotempoleOutputEntry{
	typedef TotempoleOutputEntry self_type;

public:
	TotempoleOutputEntry() :
		byte_offset(0),
		entries(0),
		uncompressed_size(0),
		contigID(-1),
		minPosition(0),
		maxPosition(0)
	{}
	~TotempoleOutputEntry(){}

	inline const bool isValid(void) const{ return(this->byte_offset != 0); }
	inline void operator++(void){ ++this->entries; }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.byte_offset << '\t' << entry.entries << '\t'
				<< entry.uncompressed_size << '\t' << entry.contigID << '\t'
				<< entry.minPosition << '\t' << entry.maxPosition;
		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.entries),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uncompressed_size), sizeof(U32));
		stream.write(reinterpret_cast<char*>(&entry.contigID),    sizeof(S32));
		stream.write(reinterpret_cast<char*>(&entry.minPosition), sizeof(U32));
		stream.write(reinterpret_cast<char*>(&entry.maxPosition), sizeof(U32));

		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.entries),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uncompressed_size), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.contigID),    sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.minPosition), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition), sizeof(U32));

		return(stream);
	}

	void reset(void){
		this->byte_offset = 0;
		this->entries = 0;
		this->uncompressed_size = 0;
		this->contigID = -1;
		this->minPosition = 0;
		this->maxPosition = 0;
	}

public:
	U64 byte_offset;	// tellg() position in stream for start of record in Tomahawk file
	U32 entries; 		// number of variants in this block
	U32 uncompressed_size; // uncompressed size of this block
	S32 contigID; // if contigID >= 0 then all entries belong to this contigID
	U32 minPosition; // minPosition of entries. 0 if contigID = -1
	U32 maxPosition; // maxPosition of entries, 0 if contigID = -1
};

#endif /* TOTEMPOLEOUTPUTENTRY_H_ */
