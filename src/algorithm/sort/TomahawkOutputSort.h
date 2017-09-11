#ifndef TOMAHAWKOUTPUTSORT_H_
#define TOMAHAWKOUTPUTSORT_H_

#include <queue>

#include "../../tomahawk/TomahawkOutput/TomahawkOutputReader.h"
#include "../../tomahawk/TomahawkOutput/TomahawkOutputManager.h"
#include "TomahawkOutputSortMergeQueueContainer.h"

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
	typedef Totempole::TotempoleOutputEntry totempoly_entry;
	typedef IO::TomahawkOutputWriterIndex writer_type;

public:
	TomahawkOutputSorter(){}
	~TomahawkOutputSorter(){}

	bool sort(const std::string& input, const std::string& destinationPrefix, const U64 memory_limit);
	bool sortMerge(const std::string& input, const std::string& destinationPrefix, const U32 block_size);

private:
	two_reader_type reader;
};


}
}
}

#endif /* TOMAHAWKOUTPUTSORT_H_ */
