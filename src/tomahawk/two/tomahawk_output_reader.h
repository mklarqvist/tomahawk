#ifndef TOMAHAWKOUTPUTREADER_H_
#define TOMAHAWKOUTPUTREADER_H_

#include <vector>
#include <fstream>
#include <string>
#include <stddef.h>
#include <regex>

#include "io/basic_buffer.h"
#include "io/compression/tgzf_controller.h"
#include "support/MagicConstants.h"
#include "algorithm/open_hashtable.h"
#include "support/type_definitions.h"
#include "third_party/intervalTree.h"
#include "tomahawk/output_container.h"
#include "tomahawk/output_container_reference.h"
#include "tomahawk/two/output_entry.h"
#include "output_filter.h"
#include "index/index.h"
#include "index/footer.h"
#include "index/tomahawk_header.h"

namespace tomahawk {

class TomahawkOutputReader {
private:
	typedef TomahawkOutputReader      self_type;
	typedef io::OutputEntry           entry_type;
	typedef OutputFilter              filter_type;
	typedef OutputContainer           output_container_type;
	typedef OutputContainerReference  output_container_reference_type;
	typedef totempole::HeaderContig   contig_type;
	typedef io::TGZFHeader            tgzf_header_type;
	typedef algorithm::ContigInterval interval_type;
	typedef TomahawkHeader            header_type;
	typedef Index                     index_type;
	typedef io::BasicBuffer           buffer_type;
	typedef io::TGZFController        tgzf_controller_type;
	typedef totempole::Footer         footer_type;
	typedef algorithm::IntervalTree<interval_type, U32> tree_type;
	typedef hash::HashTable<std::string, U32> hash_table;

public:
	TomahawkOutputReader();
	~TomahawkOutputReader();

	// Accessors
	inline footer_type& getFooter(void){ return(this->footer_); }
	inline const footer_type& getFooter(void) const{ return(this->footer_); }
	inline const index_type& getIndex(void) const{ return(*this->index_); }
	inline index_type& getIndex(void){ return(*this->index_); }
	inline const header_type& getHeader(void) const{ return(this->header_); }
	inline header_type& getHeader(void){ return(this->header_); }
	inline index_type* getIndexPointer(void){ return(this->index_); }

	/**<
	 * Primary function to open a TWO file
	 * @param input Input string file location
	 * @return      Returns TRUE if basic parsing is successful or FALSE otherwise
	 */
	bool open(const std::string input);

	/**<
	 * Adds interval regions in unparsed string format.
	 * @param positions Input vector of unparsed intervals
	 * @return          Returns TRUE if parsing is successful or FALSE otherwise
	 */
	bool addRegions(std::vector<std::string>& positions);

	// Streaming functions
	/**<
	 * Seek to block at a given position and load that
	 * data into memory
	 * @param position Target block position
	 * @return         Returns TRUE upon success or FALSE otherwise
	 */
	bool seekBlock(const U32 position);

	/**<
	 * Used in parallel programming:
	 * Takes an input file stream and seeks to a given position and load that
	 * data into memory without modifying the host container
	 * @param stream   Input file stream
	 * @param position Target block position
	 * @return         Returns TRUE upon success or FALSE otherwise
	 */
	bool seekBlock(std::ifstream& stream, const U32 position) const;

	/**<
	 * Parses TWO data that has been loaded into memory after invoking
	 * either getBlock functions. This function also increments the internal
	 * position of the file handler.
	 * @param clear     Boolean set to TRUE if raw data should be cleared after invoking this function
	 * @param clear_raw Boolean set to TRUE if compressed raw data should be cleared after invoking this function
	 * @return          Returns TRUE upon success or FALSE otherwis
	 */
	int parseBlock(const bool clear = true, const bool clear_raw = true);

	/**<
	 * Used in parallel programming:
	 * Takes an input file stream and the necessary buffers and inflates `two` data without
	 * modifying the host container
	 * @param stream              Input file stream
	 * @param inflate_buffer      Support buffer for loading compressed data
	 * @param data_buffer         Output buffer for inflated `two` entries
	 * @param compression_manager Compression manager
	 * @param clear               Boolean set to TRUE if raw data should be cleared after invoking this function
	 * @return                    Returns TRUE upon success or FALSE otherwise
	 */
	int parseBlock(std::ifstream& stream, buffer_type& inflate_buffer, buffer_type& data_buffer, tgzf_controller_type& compression_manager, const bool clear = true) const;

	bool printHeader(std::ostream& stream) const;
	bool printHeader(std::ostream& stream, std::vector<std::string>& extra);

	// Access: no random access. All these functions
	//         assumes that data is loaded linearly from disk
	inline output_container_type getContainer(void){ return(output_container_type(this->data_)); }
	output_container_type getContainerVariants(const U64 n_variants);
	output_container_type getContainerBytes(const size_t l_data);
	output_container_type getContainerBlocks(const U32 n_blocks);
	inline output_container_reference_type getContainerReference(void){ return(output_container_reference_type(this->data_)); }

	// Access: requires complete (if n = 1) or partial (if n > 1) random access
	output_container_reference_type getContainerReferenceBlock(const U32 blockID);
	output_container_reference_type getContainerReferenceBlock(std::vector<U32> blocks);
	output_container_type getContainerBlock(const U32 blockID);
	output_container_type getContainerBlock(std::vector<U32> blocks);

	inline void setShowHeader(const bool yes){ this->showHeader_ = yes; }
	inline const bool getShowHeader(void) const{ return(this->showHeader_); }
	inline const bool isSorted(void) const{ return(this->index_->getController().isSorted == true); }

	// Basic operations
	bool view(void);
	bool view(const interval_type& interval);
	bool view(const std::vector<interval_type>& intervals);

	// Concatenate
	bool concat(const std::string& file_list, const std::string& output);
	bool concat(const std::vector<std::string>& files, const std::string& output);

	inline filter_type& getFilter(void){ return this->filters_; }

	bool statistics(void);
	bool aggregate(const U32 scene_x_dimension, const U32 scene_y_dimension);

private:
	bool ParseHeader(void);
	bool ParseHeaderExtend(void);

	bool __viewOnly(void);
	bool __viewFilter(void);
	bool __viewRegion(void);

	bool __checkRegionUnsorted(const entry_type& entry);
	bool __checkRegionSorted(const entry_type& entry);
	bool __concat(const std::vector<std::string>& files, const std::string& output);

	bool __addRegions(std::vector<std::string>& positions);
	//bool __ParseRegion(const std::string& region, interval_type& interval);
	bool __ParseRegion(const std::string& region, interval_type& interval) const;
	//bool __ParseRegionIndexedBlocks(void);

public:
	U64            filesize_;  // filesize
	U64            offset_end_of_data_;
	bool           showHeader_; // flag to output header or not
	bool           output_json_;
	std::ifstream  stream_;    // reader stream

	header_type    header_;
	footer_type    footer_;
	index_type*    index_;

	buffer_type          buffer_;          // input buffer
	buffer_type          data_;            // inflate buffer
	buffer_type          outputBuffer_;    // output buffer
	tgzf_controller_type tgzf_controller_; // compression controller

	filter_type filters_;	// filter parameters

	tree_type** interval_tree; // actual interval trees
	std::vector<interval_type>* interval_tree_entries; // entries for interval trees
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTREADER_H_ */
