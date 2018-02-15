#ifndef TOMAHAWKOUTPUTSORT_H_
#define TOMAHAWKOUTPUTSORT_H_

#include <thread>
#include <queue>

#include "../../io/compression/TGZFEntryIterator.h"
#include "../../tomahawk/two/TomahawkOutputReader.h"
#include "output_sort_merge_queue.h"
#include "output_sort_slave.h"

namespace Tomahawk{
namespace Algorithm{

/**<
 * Primary class for sorting `TWO` data
 */
class OutputSorter{
	typedef IO::OutputEntry                   entry_type;
	typedef TomahawkOutputReader              two_reader_type;
	typedef IO::TGZFEntryIterator<entry_type> tgzf_iterator;
	typedef OutputSorter                      self_type;
	typedef OutputSortMergeQueue<entry_type>  queue_entry;
	typedef std::priority_queue<queue_entry>  queue_type; // priority queue

public:
	OutputSorter() : n_threads(std::thread::hardware_concurrency()){}
	~OutputSorter(){}

	/**<
	 * Standard sorting approach
	 * @param input
	 * @param destinationPrefix
	 * @param memory_limit
	 * @return
	 */
	bool sort(const std::string& input, const std::string& destinationPrefix, U64 memory_limit);

	/**<
	 * N-way merge of paralell-sorted blocks
	 * @param input
	 * @param destinationPrefix
	 * @param block_size
	 * @return
	 */
	bool sortMerge(const std::string& input, const std::string& destinationPrefix, const U32 block_size);

	inline const size_t size(void) const{ return(this->n_threads); }

private:
	two_reader_type reader;

public:
	size_t n_threads;
};


}
}

#endif /* TOMAHAWKOUTPUTSORT_H_ */
