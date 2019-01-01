#ifndef TWK_INTERVALS_H_
#define TWK_INTERVALS_H_

#include "index.h"
#include "header.h"
// Move out
#include "third_party/intervalTree.h"

namespace tomahawk {

struct interval_pair_payload {
	interval_pair_payload();
	interval_pair_payload(const int32_t rid, const uint32_t mate_offset, const uint8_t type);
	bool operator<(const interval_pair_payload& other) const;

	uint8_t mate; // right or left mate if using a linked interval A:B
	int32_t rid;
	uint32_t offset;
};

/**<
 * Basic intervals container. Handles tuples of (rid,from_pos,to_pos) and maps
 * them to the local tomahawk index. Have functionality for sorting, deduping,
 * and mapping interval tuples and the resulting index objects.
 */

class twk_intervals {
private:
	typedef algorithm::Interval<uint32_t,uint32_t> interval;

public:
	twk_intervals();
	twk_intervals(const uint32_t n_contigs);
	twk_intervals(const twk_intervals& other) = delete;
	twk_intervals& operator=(const twk_intervals& other) = delete;
	~twk_intervals();

	/**<
	 * Dedupes interval vectors tuples (rid,fromA,fromB) by extension. If two
	 * intervals partially overlap in the right end then we extend the left-most
	 * interval to the right-most intervals end position.
	 */
	void Dedupe();

	/**<
	 * Wrapper function for mapping the internal interval tuples to the index
	 * entries provided by the twk_reader object. If construction fails or
	 * no map hits are found we return FALSE.
	 * @param reader Source twk_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool Build(const uint32_t n_contigs, const Index& index);

	/**<
	 * Convenince wrapper for interval strings. Iterate over a vector of unparsed
	 * interval strings and parse them.
	 * @param ivals  Source vector of unparsed interval strings.
	 * @param reader Reference instance of twk_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseIntervalStrings(std::vector<std::string>& ivals, VcfHeader& hdr);

	/**<
	 * Internal function for parsing a source interval string into a (rid,posA,posB)
	 * tuple. Input string has to match any of the following patterns:
	 *    1) ^rid$
	 *    2) ^rid:pos$
	 *    3) ^rid:pos-pos$
	 *
	 * Numerical values may be expressed in any standard form: for example, 10e6
	 * and 20E6 and 100000. This function will return TRUE if parsing is syntactically
	 * legal and the contig identifier exists in the header. This function will not check
	 * for out-of-bounds errors such as providing an interval exceeding the length of a
	 * contig.
	 * @param s      Source unparsed interval string/
	 * @param reader Reference instance of twk_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseIntervalString(const std::string& s, const VcfHeader& hdr);

public:
	uint32_t n_c;
	std::vector< std::vector< interval > > ivecs; // vector of vectors of intervals
	algorithm::IntervalTree<uint32_t,uint32_t>** itree; // interval tree array
	std::vector<IndexEntry*> overlap_blocks; // overlapping blocks of interest
};

class twk_intervals_two {
public:
	// interval (from,to,offset to mate) with rid being implicit in the tree
	typedef algorithm::Interval< uint32_t, interval_pair_payload > interval;

	twk_intervals_two();
	twk_intervals_two(const uint32_t n_contigs);
	~twk_intervals_two();

	/**<
	 * Wrapper function for mapping the internal interval tuples to the index
	 * entries provided by the twk_reader object. If construction fails or
	 * no map hits are found we return FALSE.
	 * @param reader Source twk_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool Build(std::vector<std::string>& strings,
	           const uint32_t n_contigs,
	           const IndexOutput& index,
	           const VcfHeader& hdr);

	/**<
	 * Internal function for parsing a source interval string into a (rid,posA,posB)
	 * tuple. Input string has to match any of the following patterns:
	 *    1) ^rid$
	 *    2) ^rid:pos$
	 *    3) ^rid:pos-pos$
	 *    4) ^rid;rid$
	 *    5) ^rid:pos,rid$
	 *    6) ^rid:pos-pos,rid$
	 *    7) ^rid,rid:pos$
	 *    8) ^rid,rid:pos-pos$
	 *
	 * Numerical values may be expressed in any standard form: for example, 10e6
	 * and 20E6 and 100000. This function will return TRUE if parsing is syntactically
	 * legal and the contig identifier exists in the header. This function will not check
	 * for out-of-bounds errors such as providing an interval exceeding the length of a
	 * contig.
	 * @param s      Source unparsed interval string
	 * @param reader Reference instance of two_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseIntervalString(const std::string& s, const VcfHeader& hdr);

	/**<
	 * Dedupes interval vectors tuples (rid,fromA,fromB) by extension. If two
	 * intervals partially overlap in the right end then we extend the left-most
	 * interval to the right-most intervals end position.
	 */
	void Dedupe();

	/**<
	 * Predicate for filtering out a provided twk1_two_t record.
	 * @param rec Input reference two record.
	 * @return    Returns TRUE if filtered out or FALSE otherwise.
	 */
	bool FilterInterval(const twk1_two_t& rec) const;

	// Accessors
	size_t GetOverlapSize() const{ return(this->overlap_blocks.size()); }
	IndexEntryOutput* GetOverlapBlock(const uint32_t p){ return(overlap_blocks[p]); }

private:
	uint32_t n_c;
	std::vector< std::vector< interval > > ivecs; // vector of vectors of intervals
	std::vector< std::vector< interval > > ivecs_internal; // deduped for finding overlapping indices
	algorithm::IntervalTree< uint32_t, interval_pair_payload >** itree; // interval tree array
	std::vector<IndexEntryOutput*> overlap_blocks; // overlapping blocks of interest
};

}


#endif /* TWK_INTERVALS_H_ */
