#ifndef TOMAHAWKOUTPUTSORT_H_
#define TOMAHAWKOUTPUTSORT_H_

#include <queue>

#include "../../tomahawk/TomahawkOutput/TomahawkOutputReader.h"
#include "TomahawkOutputSortSupport.h"
#include "../../tomahawk/TomahawkOutput/TomahawkOutputManager.h"

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

// Sorter
class TomahawkOutputSorter{
	typedef IO::TomahawkOutputEntry entry_type;
	typedef IO::TomahawkOutputEntrySort entry_sort_type;
	typedef TomahawkOutputSorter self_type;
	typedef TomahawkOutputSortMergeQueueContainer<entry_type> queue_entry;
	typedef std::priority_queue< queue_entry > queue_type; // prio queue
	typedef IO::TomahawkOutputReader two_reader_type;
	typedef IO::PartialSortIndexHeader partial_header_type;
	typedef IO::PartialSortIndexHeaderEntry partial_header_entry_type;

public:
	TomahawkOutputSorter(){}
	~TomahawkOutputSorter(){}

	//bool sort(const std::string& input, const U64 memory_limit);
	bool sort(const std::string& input, const std::string& destinationPrefix, const U64 memory_limit);
	bool sortMerge(const std::string& input);

private:
	template <class S> bool sortMerge(const U32 count, S& outstream);

private:
	two_reader_type reader;
};


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
