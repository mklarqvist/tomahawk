#ifndef PACKEDENTRYREADER_H_
#define PACKEDENTRYREADER_H_

namespace Tomahawk{
namespace IO{

#define PACKED_READER_DEFAULT_CHUNK	1000000 // 1MB

template <class T, int Y = sizeof(T)>
class PackedEntryReader{
	typedef TomahawkOutputEntry entry_type;

public:
	PackedEntryReader();
	virtual ~PackedEntryReader();

	bool setup(const std::string file, size_t chunk_size = PACKED_READER_DEFAULT_CHUNK - PACKED_READER_DEFAULT_CHUNK % Y);
	bool nextEntry(const entry_type*& entry);
	virtual bool nextBlock(void);

	inline bool seek(const size_t pos);
	inline entry_type* begin(void){ return(this->entries); }
	inline entry_type* end(void){ return(&this->entries[this->entry_tail-1]); }
	inline entry_type* operator[](const U32& p){ return(&this->entries[p]); }
	inline const size_t& size(void) const{ return this->entry_tail; }
	inline const size_t& size_buffer(void) const{ return this->buffer_size; }
	inline void reset(void){ this->entry_head = 0; this->entry_tail = 0; }
	inline void next(void){ ++this->entry_head; }
	inline void prev(void){ --this->entry_head; }
	inline bool available(void) const{ return(this->entry_head < this->entry_tail); }
	inline bool good(void) const{ return(this->stream.good()); }
	inline const size_t& filesize(void) const{ return this->__filesize; }
	inline size_t tellg(void){ return this->stream.tellg(); }
	inline const size_t& block_size(void) const{ return this->read_block_size; }

protected:
	bool open(const std::string& file);

protected:
	size_t __filesize;
	size_t entry_head;
	size_t entry_tail;
	size_t buffer_size;
	size_t read_block_size;
	std::ifstream stream;
	char* buffer;
	entry_type* entries;
};

template <class T, int Y>
PackedEntryReader<T, Y>::PackedEntryReader()
	: __filesize(0)
	, entry_head(0)
	, entry_tail(0)
	, buffer_size(0)
	, read_block_size(0)
	, buffer(nullptr)
	, entries(nullptr)
{}

template <class T, int Y>
PackedEntryReader<T, Y>::~PackedEntryReader(){
	delete [] this->buffer;
}

template <class T, int Y>
bool PackedEntryReader<T, Y>::setup(const std::string file, size_t chunk_size){
	if(!this->open(file)){
		std::cerr << Helpers::timestamp("ERROR", "IO") << "Failed to open file..." << std::endl;
		return false;
	}

	if(chunk_size == 0){
		std::cerr << "illegal chunk size" << std::endl;
		return false;
	}

	if(chunk_size % Y != 0){
		std::cerr << "Adjusting chunk size: " << chunk_size << " -> ";
		chunk_size -= chunk_size % Y;
		std::cerr << chunk_size << std::endl;
		if(chunk_size == 0){
			std::cerr << "illegal chunk size" << std::endl;
			return false;
		}
	}

	this->read_block_size = chunk_size;

	this->reset();
	delete [] this->buffer;
	this->buffer = new char[this->read_block_size];
	this->buffer_size = this->read_block_size;

	return true;
}

template <class T, int Y>
bool PackedEntryReader<T, Y>::open(const std::string& file){
	this->stream.open(file, std::ios::binary | std::ios::in | std::ios::ate);
	if(!this->good()){
		std::cerr << Helpers::timestamp("ERROR", "IO") << "IO-stream is bad..." << std::endl;
		return false;
	}

	this->__filesize = this->stream.tellg();
	this->stream.seekg(0);

	return true;
}

template <class T, int Y>
bool PackedEntryReader<T, Y>::nextEntry(const entry_type*& entry){
	if(!this->available()){
		if(!this->nextBlock())
			return false;

	}

	entry = &this->entries[this->entry_head];
	this->next();
	return true;
}

template <class T, int Y>
bool PackedEntryReader<T, Y>::nextBlock(void){
	if(!this->good()){
		std::cerr << Helpers::timestamp("ERROR", "IO") << "IO-stream has failed..." << std::endl;
		return false;
	}

	if(this->stream.tellg() == this->filesize())
		return false;

	// Ignore if unset
	// comparison does not happen if tellg() == -1, return above
	this->reset();

	size_t readAmount = this->read_block_size;
	if((U64)this->stream.tellg() + this->read_block_size > this->filesize())
		readAmount = this->filesize() - this->stream.tellg();

	this->stream.read(this->buffer, readAmount);
	const U32 entries_read = this->stream.gcount() / Y;
	this->buffer_size = this->stream.gcount();
	if(this->stream.gcount() % Y != 0){
		std::cerr << Helpers::timestamp("ERROR", "IO") << "block is staggered" << std::endl;
		return false;
	}

	this->entry_tail = entries_read;
	this->entries = reinterpret_cast<entry_type*>(this->buffer);

	return true;
}

template <class T, int Y>
bool PackedEntryReader<T, Y>::seek(const size_t pos){
	this->stream.seekg(pos);
	this->reset(); // trigger reloading data when asking for next entry
	return(this->good());
}

}
}



#endif /* PACKEDENTRYREADER_H_ */
