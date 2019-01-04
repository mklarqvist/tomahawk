/*
Copyright (C) 2016-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TWK_LD_H_
#define TWK_LD_H_

#include <cassert>
#include <thread>

#include "core.h"
#include "twk_reader.h"
#include "writer.h"
// Move out
#include "intervals.h"

namespace tomahawk {

struct twk_ld_balancer; // forward declare balancer

/****************************
*  LD handler
****************************/
class twk_ld {
public:
	twk_ld();
	~twk_ld();

	/**<
	 * Wrapper function for parsing interval strings.
	 * @param reader Reference twk_reader instance.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseIntervalStrings(twk_reader& reader){ return(intervals.ParseIntervalStrings(settings.ival_strings, reader.hdr)); }

	void operator=(const twk_ld_settings& settings){ this->settings = settings; }

	/**<
	 * Reads the desired tomahawk blocks given the balancer intervals. Internally
	 * decides if the block slicing is based on the universal set of blocks or a
	 * targetted subset (as decided by the interval slicing operation).
	 * @param reader
	 * @param bit
	 * @param balancer
	 * @param load
	 * @return
	 */
	bool LoadBlocks(twk_reader& reader,
	                twk1_blk_iterator& bit,
	                const twk_ld_balancer& balancer,
	                const twk_ld_settings& settings);

	/**<
	 * Loading twk blocks for a single variant and its surrounding variants within
	 * some distance and on the same chromosome. This function loads the target variant
	 * in a single block and the other variants in separate blocks as usual. The
	 * identity of the target site is parameterized in the settings object.
	 * @param reader
	 * @param bit
	 * @param balancer
	 * @param load
	 * @return
	 */
	bool LoadTargetSingle(twk_reader& reader,
	                      twk1_blk_iterator& bit,
	                      const twk_ld_balancer& balancer,
	                      const uint8_t load);

	/**<
	 * Construct interval container and trees given the pre-provided interval
	 * strings.
	 * @param reader
	 * @param bit
	 * @param load
	 * @return
	 */
	bool BuildIntervals(twk_reader& reader, twk1_blk_iterator& bit);

	/**<
	 * Loads only the target blocks that overlap with the given vector of interval
	 * tuples as parameterized in the settings object.
	 * @param reader   Reference to twk reader.
	 * @param bit      Reference to a twk block iterator.
	 * @param balancer Reference of a pre-computed load balancer.
	 * @param settings Reference of a user-paramterized settings object.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool LoadTargetBlocks(twk_reader& reader,
	                      twk1_blk_iterator& bit,
	                      const twk_ld_balancer& balancer,
	                      const twk_ld_settings& settings);

	/**<
	 * Loads all available twk blocks into memory. Internally spawns the maximum
	 * possible number of unpacking threads possible (as parameterized by settings).
	 * @param reader   Reference to twk reader.
	 * @param bit      Reference to a twk block iterator.
	 * @param balancer Reference of a pre-computed load balancer.
	 * @param settings Reference of a user-paramterized settings object.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool LoadAllBlocks(twk_reader& reader,
	                   twk1_blk_iterator& bit,
	                   const twk_ld_balancer& balancer,
	                   const twk_ld_settings& settings);

	/**<
	 * Helper function to call Compute subroutine when passing a new settings
	 * object to be used.
	 * @param settings Src settings objects.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool Compute(const twk_ld_settings& settings);

	/**<
	 * Main subroutine for computing linkage-disequilibrium as contextually
	 * determined given the user-defined parameters.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool Compute();
	bool ComputePerformance();

private:
	class twk_ld_impl;
	uint32_t n_blks, m_blks, n_vnts, n_tree;
	twk1_block_t* ldd2;
	twk_ld_settings settings;
	twk_intervals intervals;
	twk_ld_impl* mImpl;
};

}

#endif
