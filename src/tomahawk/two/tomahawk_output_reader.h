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
#include "tomahawk/two/aggregation_parameters.h"
#include "output_filter.h"
#include "index/index.h"
#include "index/footer.h"
#include "index/tomahawk_header.h"
#include "io/output_writer.h"

namespace tomahawk {

enum TOMAHAWK_OUTPUT_FAMILY{
	TWK_OUTPUT_TWO,
	TWK_OUTPUT_LD
};

struct TomahawkOutputReaderParameters{
public:
	TomahawkOutputReaderParameters() :
		showHeader(true),
		output_json(false),
		output_type(TWK_OUTPUT_LD),
		output_file("-")
	{

	}

	~TomahawkOutputReaderParameters() = default;

public:
	bool showHeader;
	bool output_json;
	TOMAHAWK_OUTPUT_FAMILY output_type;
	std::string input_file;
	std::string output_file;
};

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
	typedef TomahawkOutputReaderParameters parameters_type;

	typedef io::OutputWriterInterface    writer_type;
	typedef io::OutputWriterBinaryStream writer_binary_stream_type;
	typedef io::OutputWriterBinaryFile   writer_binary_file_type;
	typedef io::OutputWriterStdOut       writer_ld_stream_type;

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

	inline const bool isSorted(void) const{ return(this->index_->getController().isSorted == true); }

	// Basic operations
	bool view(void);
	//bool view(const interval_type& interval);
	//bool view(const std::vector<interval_type>& intervals);

	// Concatenate
	bool concat(const std::string& file_list, const std::string& output);
	bool concat(const std::vector<std::string>& files, const std::string& output);

	inline filter_type& getFilter(void){ return this->filters_; }

	bool statistics(void);
	bool aggregate(support::aggregation_parameters& parameters);

private:
	inline bool openWriter(void){
		if(this->parameters_.output_type == TWK_OUTPUT_TWO){
			if(this->parameters_.output_file == "-" || this->parameters_.output_file.size() == 0){
				this->writer_ = new writer_binary_stream_type;
				return(true);
			} else {
				this->writer_ = new writer_binary_file_type;
				if(!this->writer_->open(this->parameters_.output_file)){
					std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to open output file handle: " << this->parameters_.output_file << std::endl;
					return false;
				}
				return(true);
			}
		} else {
			this->writer_ = new writer_ld_stream_type;
			return(true);
		}
		return(false);
	}

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
	std::ifstream  stream_;    // reader stream

	parameters_type parameters_;
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
	writer_type* writer_;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTREADER_H_ */
