#ifndef TOMAHAWKOUTPUTSORT_H_
#define TOMAHAWKOUTPUTSORT_H_

#include <queue>

#include "TomahawkOutputSortSupport.h"

namespace Tomahawk{
namespace Algorithm{
namespace Output{

template <class T>
struct TomahawkOutputSortMergeQueueContainer {
	typedef T entry_type;
	typedef TomahawkOutputSortMergeQueueContainer<entry_type> self_type;

public:
	TomahawkOutputSortMergeQueueContainer(const entry_type* data,
										  U32 streamID,
										  bool (*compFunc)(const entry_type* a, const entry_type* b) = T::operator<)
	: streamID(streamID)
	, data(data)
	, compFunc(compFunc)
    {}

    bool operator<(const self_type& a) const{
        return!(this->compFunc(this->data, a.data));
    }

public:
    U32 streamID;
    const entry_type* data;
    bool (*compFunc)(const entry_type* a, const entry_type* b);
};

// This class extends the PackedEntryReader class
// as it requires minor changes to how setup
// and retrieving next block is handled
// because of IO bounds (start and end positions)
template <class T, int Y = sizeof(T)>
class TomahawkOutputSortReader : public Tomahawk::IO::PackedEntryReader<T, Y> {
public:
	typedef T entry_type;

public:
	TomahawkOutputSortReader();
	~TomahawkOutputSortReader();

	bool setup(const std::string file, const U64 start_offset, const U64 upper_bounds, const U64 chunk_size);
	bool nextBlock(void);

	inline U64& getBounds(void) const{ return this->stream_bounds; }

private:
	bool open(const std::string& file, const U64 start_offset);

private:
	U64 stream_bounds; // where stream should terminate
};

template <class T, int Y>
TomahawkOutputSortReader<T, Y>::TomahawkOutputSortReader(void)
	: stream_bounds(0)
{}

template <class T, int Y>
TomahawkOutputSortReader<T, Y>::~TomahawkOutputSortReader(void){}

template <class T, int Y>
bool TomahawkOutputSortReader<T, Y>::setup(const std::string file, const U64 start_offset, const U64 upper_bounds, const U64 chunk_size){
	if(!this->open(file, start_offset)){
		std::cerr << "failed to open " << file << " at " << start_offset << std::endl;
		return false;
	}

	this->stream_bounds = upper_bounds;
	this->read_block_size = chunk_size;

	this->reset();
	delete [] this->buffer;
	this->buffer = new char[this->read_block_size];

	return true;
}

template <class T, int Y>
bool TomahawkOutputSortReader<T, Y>::open(const std::string& file, const U64 start_offset){
	this->stream.open(file, std::ios::binary | std::ios::in);
	if(!this->good()) return false;
	this->stream.seekg(start_offset);
	if(!this->good()) return false;

	return true;
}


template <class T, int Y>
bool TomahawkOutputSortReader<T, Y>::nextBlock(void){
	if(!this->good()){
		std::cerr << "no good" << std::endl;
		return false;
	}

	// Ignore if unset
	// comparison does not happen if tellg() == -1, return above
	if(this->stream.tellg() == this->stream_bounds){
		//std::cerr << "at bounds: " << this->stream.tellg() << '/' << this->stream_bounds << std::endl;
		return false;
	}

	if(this->stream.tellg() >= this->stream_bounds){
		std::cerr << "out of bounds: " << this->stream.tellg() << '/' << this->stream_bounds << std::endl;
		return false;
	}

	this->reset();
	const U64 testBounds = (U64)this->stream.tellg() + this->read_block_size;
	size_t readBlockSize = this->read_block_size;
	if(testBounds > this->stream_bounds) readBlockSize = this->stream_bounds - this->stream.tellg();

	this->stream.read(this->buffer, readBlockSize);
	const U32 entries_read = this->stream.gcount() / Y;
	if(this->stream.gcount() % Y != 0){
		std::cerr << "block is staggered" << std::endl;
		return false;
	}

	this->entry_tail = entries_read;
	this->entries = reinterpret_cast<entry_type*>(this->buffer);

	return true;
}

// Sorter
template <class T, int Y = sizeof(T)>
class TomahawkOutputSorter{
	typedef T entry_type;
	typedef TomahawkOutputSorter<entry_type, Y> self_type;
	typedef IO::PackedEntryReader<entry_type, Y> packed_reader;
	typedef TomahawkOutputSortReader<entry_type, Y> sort_reader;
	typedef TomahawkOutputSortMergeQueueContainer<entry_type> queue_entry;
	typedef std::priority_queue< queue_entry > queue_type;

	typedef IO::PartialSortIndexHeader partial_header_type;
	typedef IO::PartialSortIndexHeaderEntry partial_header_entry_type;

public:
	TomahawkOutputSorter();
	~TomahawkOutputSorter();

	bool sort(const std::string& input);
	bool sort(const std::string& input, const std::string& destinationPrefix);
	bool kwayMerge(const std::string& input);

private:
	bool sort(const std::string& input, const std::string& outFile, const std::string& indexOut);

	// Type S requires function write(char*, length)
	template <class S> bool sort(packed_reader& reader, S& outstream, S& indexstream);
	template <class S> bool sortInplace(S& outstream);
	template <class S> bool kwayMerge(sort_reader* reader, const U32 count, S& outstream);
};

#pragma pack(1)
struct TomahawkOutputSortIndexEntry{
	typedef TomahawkOutputSortIndexEntry self_type;

public:
	friend std::ostream& operator<<(std::ostream& os, const self_type& entry){
		os << entry.from << '\t' << entry.to;
		return(os);
	}

	friend std::ofstream& operator<<(std::ofstream& of, const self_type& entry){
		of.write((char*)&entry.from, sizeof(U64));
		of.write((char*)&entry.to,   sizeof(U64));
		return(of);
	}

public:
	U64 from;
	U64 to;
};

template <class T, int Y>
TomahawkOutputSorter<T,Y>::TomahawkOutputSorter(){}

template <class T, int Y>
TomahawkOutputSorter<T,Y>::~TomahawkOutputSorter(){}

template <class T, int Y>
bool TomahawkOutputSorter<T,Y>::sort(const std::string& input){
	std::vector<std::string> paths = Tomahawk::Helpers::splitLastOf(input, '/', true);
	std::vector<std::string> files = Tomahawk::Helpers::splitLastOf(paths[1], '.');

	// Todo: if failed to read from file suffix: try to look into file header MAGIC
	if(files[1].size() == 0){
		std::cerr << "could not determine file type from suffix" << std::endl;
		return false;
	}

	std::transform(files[1].begin(), files[1].end(), files[1].begin(), ::tolower);

	// Setup filenames
	std::string outFile = files[0] + "__partial_sort." + Tomahawk::Constants::OUTPUT_LD_SUFFIX;
	std::string outIndex = outFile + '.' + Tomahawk::Constants::OUTPUT_LD_PARTIAL_SORT_INDEX_SUFFIX;
	std::cerr << Helpers::timestamp("LOG","SORT") << "Setting filename to: " << outFile << std::endl;
	std::cerr << Helpers::timestamp("LOG","SORT") << "Setting index filename to: " << outIndex << std::endl;
	outFile = paths[0] + outFile;
	outIndex = paths[0] + outIndex;

	return(this->sort(input, outFile, outIndex));
}

template <class T, int Y>
bool TomahawkOutputSorter<T,Y>::sort(const std::string& input, const std::string& destinationPrefix){
	std::string outFile = destinationPrefix + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX;
	std::string outIndex = outFile + '.' + Tomahawk::Constants::OUTPUT_LD_PARTIAL_SORT_INDEX_SUFFIX;
	return(this->sort(input, outFile, outIndex));
}

template <class T, int Y>
bool TomahawkOutputSorter<T,Y>::sort(const std::string& input, const std::string& outFile, const std::string& indexOut){
	U64 blockSize = 1e9; // 1GB blocks
	blockSize -= blockSize % Y;

	IO::PackedEntryReader<T, Y> reader;
	if(!reader.setup(input, blockSize))
		return false;

	std::ofstream outFileStream(outFile, std::ios::out | std::ios::binary);
	if(!outFileStream.good()){
		std::cerr << "Faield to open outfile" << std::endl;
		return false;
	}

	std::ofstream outIndexStream(indexOut, std::ios::out | std::ios::binary);
	if(!outIndexStream.good()){
		std::cerr << "Faield to open outindex" << std::endl;
		return false;
	}
	outIndexStream.write(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH);
	outIndexStream.write((char*)&Tomahawk::Constants::PROGRAM_VERSION, sizeof(float));

	return(this->sort(reader, outFileStream, outIndexStream));
}

template <class T, int Y>
template <class S>
bool TomahawkOutputSorter<T,Y>::sort(packed_reader& reader, S& outstream, S& indexstream){
	const U32 entries = reader.block_size()/Y;
	std::cerr << Helpers::timestamp("LOG","SORT") << "Sorting " << reader.filesize() << " in blocks of " << reader.block_size() << "..." << std::endl;

	// Pointers to reinterpreted pointers used for sorting
	entry_type** pointers = new entry_type*[entries];

	// Output partial sort index
	TomahawkOutputSortIndexEntry index_entry;
	index_entry.from = 0;

	while(reader.nextBlock()){
		index_entry.to = reader.tellg();
		indexstream << index_entry;
		std::swap(index_entry.from, index_entry.to);

		for(U32 i = 0; i < reader.size(); ++i)
			pointers[i] = reader[i];

		std::sort(pointers, &pointers[reader.size()], Tomahawk::IO::Support::TomahawkOutputEntryCompFunc);
		for(U32 i = 0; i < reader.size(); ++i)
			outstream.write(reinterpret_cast<const char*>(pointers[i]), sizeof(entry_type));
	}

	// Flush output (should be done automatically upon leaving function)
	outstream.flush();
	indexstream.flush();

	delete [] pointers;

	return true;
}

template <class T, int Y>
bool TomahawkOutputSorter<T,Y>::kwayMerge(const std::string& inputFile){
	std::cerr << "attempting to open: " << inputFile + '.' + Tomahawk::Constants::OUTPUT_LD_PARTIAL_SORT_INDEX_SUFFIX << std::endl;
	std::ifstream indexStream(inputFile + '.' + Tomahawk::Constants::OUTPUT_LD_PARTIAL_SORT_INDEX_SUFFIX, std::ios::binary | std::ios::in | std::ios::ate);
	if(!indexStream.good()){
		std::cerr << "bad index stream" << std::endl;
		return false;
	}
	const U64 filesize_index = indexStream.tellg();
	indexStream.seekg(0);

	// parse header
	partial_header_type head;
	indexStream >> head;
	if(!head.validate()){
		std::cerr << "Failed to validate header" << std::endl;
		return false;
	}
	std::cerr << std::string(&head.header[0], Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH) << '\t' << head.version << std::endl;
	std::cerr << "position now: " << indexStream.tellg() << '/' << filesize_index << std::endl;
	std::cerr << "entries: " << (filesize_index - indexStream.tellg())/(2*sizeof(U64)) << std::endl;

	char* indexHeaderEntries = new char[filesize_index - indexStream.tellg()];
	const U32 indexEntryEnd = (filesize_index - indexStream.tellg())/(2*sizeof(U64));
	indexStream.read(indexHeaderEntries, filesize_index - indexStream.tellg());
	const partial_header_entry_type* const entries = reinterpret_cast<const partial_header_entry_type* const>(&indexHeaderEntries[0]);

	sort_reader* sortEntries = new sort_reader[indexEntryEnd];

	std::cerr << Helpers::timestamp("LOG", "SORT") << "Opening " << indexEntryEnd << " file handles..." << std::endl;
	for(U32 i = 0; i < indexEntryEnd; ++i){
		//std::cerr << entries[i] << std::endl;

		if(!sortEntries[i].setup(inputFile, entries[i].from, entries[i].to, 1000000 - (1000000 % sizeof(entry_type)))){
			std::cerr << "failed setup" << std::endl;
		}
	}

	//
	if(!this->kwayMerge(sortEntries, indexEntryEnd, std::cout)){
		delete [] sortEntries;
		delete [] indexHeaderEntries;
		return false;
	}

	delete [] sortEntries;
	delete [] indexHeaderEntries;

	return true;
}

template <class T, int Y>
template <class S>
bool TomahawkOutputSorter<T,Y>::kwayMerge(sort_reader* reader, const U32 count, S& outstream){
	// queue
	queue_type outQueue;

	// draw one from each
	const entry_type* e;
	for(U32 i = 0; i < count; ++i){
		reader[i].nextEntry(e);
		outQueue.push( queue_entry(e, i, IO::Support::TomahawkOutputEntryCompFuncConst) );
	}

	 // while queue is not empty
	while(outQueue.empty() == false){
		// peek at top entry in queue
		const queue_entry& lowest = outQueue.top();
		const U32 id = lowest.streamID;
		// write the entry from the top of the queue
		std::cout.write(reinterpret_cast<const char*>(lowest.data), Y);
		// remove this record from the queue
		outQueue.pop();

		// If it is possible to get a new entry from this particular stream
		if(reader[id].nextEntry(e)){
			// Push new entry into priority queue
			outQueue.push( queue_entry(e, id, IO::Support::TomahawkOutputEntryCompFuncConst) );
		}
	}

	return true;
}


}
}
}



#endif /* TOMAHAWKOUTPUTSORT_H_ */
