#ifndef TOMAHAWKOUTPUTSORT_H_
#define TOMAHAWKOUTPUTSORT_H_

#include <queue>

#include "../../io/TomahawkOutput/TomahawkOutputReader.h"
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

// Sorter
class TomahawkOutputSorter{
	typedef IO::TomahawkOutputEntry entry_type;
	typedef TomahawkOutputSorter self_type;
	typedef TomahawkOutputSortMergeQueueContainer<entry_type> queue_entry;
	typedef std::priority_queue< queue_entry > queue_type; // prio queue
	typedef IO::TomahawkOutputReader two_reader_type;
	typedef IO::PartialSortIndexHeader partial_header_type;
	typedef IO::PartialSortIndexHeaderEntry partial_header_entry_type;

public:
	TomahawkOutputSorter(){}
	~TomahawkOutputSorter(){}

	bool sort(const std::string& input);
	bool sort(const std::string& input, const std::string& destinationPrefix);
	bool sortMerge(const std::string& input);

private:
	bool sort(const std::string& input, const std::string& outFile, const std::string& indexOut);

	// Type S requires function write(char*, length)
	template <class S> bool sort(S& outstream, S& indexstream);
	//template <class S> bool sortInplace(S& outstream);
	template <class S> bool sortMerge(const U32 count, S& outstream);

private:
	two_reader_type reader;

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



template <class S>
bool TomahawkOutputSorter::sort(S& outstream, S& indexstream){
	/*
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
	*/
}

template <class S>
bool TomahawkOutputSorter::sortMerge(const U32 count, S& outstream){
	/*
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
	*/
}


}
}
}

#endif /* TOMAHAWKOUTPUTSORT_H_ */
