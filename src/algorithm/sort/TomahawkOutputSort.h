#ifndef TOMAHAWKOUTPUTSORT_H_
#define TOMAHAWKOUTPUTSORT_H_

#include <thread>
#include <queue>

#include "../../io/compression/TGZFEntryIterator.h"
#include "../../totempole/TotempoleOutputEntry.h"
#include "../../tomahawk/TomahawkOutput/TomahawkOutputReader.h"
#include "../../tomahawk/TomahawkOutput/TomahawkOutputManager.h"
#include "TomahawkOutputSortMergeQueueContainer.h"
#include "TomahawkOutputSortSlave.h"

namespace Tomahawk{
namespace Algorithm{
namespace Output{

// Sorter
class TomahawkOutputSorter{
	typedef IO::TomahawkOutputEntry entry_type;
	typedef IO::TomahawkOutputEntrySort entry_sort_type;
	typedef TomahawkOutputSorter self_type;
	typedef TomahawkOutputSortMergeQueueContainer<entry_type> queue_entry;
	typedef std::priority_queue< queue_entry > queue_type; // prio queue
	typedef IO::TomahawkOutputReader two_reader_type;
	typedef Totempole::TotempoleOutputEntry totempole_entry;
	typedef IO::TomahawkOutputWriterIndex writer_type;
	typedef IO::WriterFile basic_writer_type;
	typedef TomahawkOutputSortSlave slave_sorter;
	typedef IO::TGZFEntryIterator<entry_type> tgzf_iterator;
	typedef Totempole::TotempoleOutputSortedIndex index_type;
	typedef Tomahawk::IO::TomahawkOutputSortHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> toi_header_type;

public:
	TomahawkOutputSorter() : n_threads(std::thread::hardware_concurrency()), reverse_entries(true){}
	~TomahawkOutputSorter(){}

	bool sort(const std::string& input, const std::string& destinationPrefix, U64 memory_limit);
	bool sortMerge(const std::string& input, const std::string& destinationPrefix, const U32 block_size);

private:
	bool __sortUnindexed();
	bool __sortIndexed(basic_writer_type& toi_writer, const std::string& input, U64 memory_limit);

private:
	two_reader_type reader;

public:
	U32 n_threads;
	bool reverse_entries;
	std::string baseName;
	std::string basePath;
};


}
}
}

#endif /* TOMAHAWKOUTPUTSORT_H_ */
