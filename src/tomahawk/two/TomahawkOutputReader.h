#ifndef TOMAHAWKOUTPUTREADER_H_
#define TOMAHAWKOUTPUTREADER_H_

#include <vector>
#include <fstream>
#include <string>
#include <stddef.h>
#include <regex>

#include "../../io/BasicBuffer.h"
#include "../../io/compression/TGZFController.h"
#include "../../support/MagicConstants.h"
#include "../../algorithm/OpenHashTable.h"
#include "../../support/type_definitions.h"
#include "../../totempole/TotempoleMagic.h"
#include "../../third_party/intervalTree.h"
#include "../../totempole/TotempoleOutputReader.h"
#include "../../tomahawk/output_container.h"
#include "../../tomahawk/output_container_reference.h"
#include "../two/output_entry.h"
#include "output_filter.h"
#include "output_writer.h"

namespace Tomahawk {
namespace IO {

class TomahawkOutputReader {
private:
	typedef OutputEntry                    entry_type;
	typedef OutputFilter                   filter_type;
	typedef OutputWriterInterface          writer_type;
	typedef OutputWriterIndex              twoi_writer_type;
	typedef OutputContainer                output_container_type;
	typedef OutputContainerReference       output_container_reference_type;
	typedef IO::BasicBuffer                buffer_type;
	typedef Totempole::TotempoleContigBase contig_type;
	typedef Totempole::TotempoleOutputSortedEntry totempole_sorted_entry_type;
	typedef TGZFHeader                     tgzf_header_type;
	typedef TGZFController                 tgzf_controller_type;
	typedef Hash::HashTable<std::string, U32> hash_table;
	typedef Algorithm::ContigInterval      interval_type;
	typedef Algorithm::IntervalTree<interval_type, U32> tree_type;
	typedef TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> header_type;
	typedef TomahawkOutputSortHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> toi_header_type;

public:
	typedef Totempole::TotempoleOutputReader toi_reader_type;
	enum WRITER_TYPE {binary, natural};

public:
	TomahawkOutputReader();
	~TomahawkOutputReader();

	const entry_type* operator[](const U32 p) const{ return(reinterpret_cast<const entry_type*>(&this->data_buffer[sizeof(entry_type)*p])); }

	bool addRegions(std::vector<std::string>& positions);
	bool Open(const std::string input);
	bool OpenExtend(const std::string input);
	inline void addLiteral(const std::string& string){ this->literals += string; }

	// Streaming functions
	/**<
	 * Seek to block at a given position and load that
	 * data into memory
	 * @param position Target block position
	 * @return         Returns TRUE upon success or FALSE otherwise
	 */
	bool seekBlock(const U32 position);

	/**<
	 * Parses TWO data that has been loaded into memory after invoking
	 * either getBlock functions. This function also increments the internal
	 * position of the file handler.
	 * @param clear Boolean set to TRUE if raw data should be cleared after invoking this function
	 * @return      Returns TRUE upon success or FALSE otherwis
	 */
	int parseBlock(const bool clear = true);

	// Access: no random access. All these functions
	//         assumes that data is loaded linearly from disk
	inline output_container_type getContainer(void){ return(output_container_type(this->data_buffer)); }
	output_container_type getContainerVariants(const U64 n_variants);
	output_container_type getContainerBytes(const size_t l_data);
	inline output_container_reference_type getContainerReference(void){ return(output_container_reference_type(this->data_buffer)); }

	// Access: requires complete (if n = 1) or partial (if n > 1) random access
	output_container_reference_type getContainerReferenceBlock(const U32 blockID);
	output_container_reference_type getContainerReferenceBlock(std::vector<U32> blocks);
	output_container_type getContainerBlock(const U32 blockID);
	output_container_type getContainerBlock(std::vector<U32> blocks);

	// Basic operations
	bool view(void);
	bool view(const interval_type& interval);
	bool view(const std::vector<interval_type>& intervals);

	// Other
	bool view(const std::string& filename);
	bool index(const std::string& filename);
	bool summary(const std::string& input, const U32 bins);

	// Concatenate
	bool concat(const std::string& file_list, const std::string& output);
	bool concat(const std::vector<std::string>& files, const std::string& output);

	//
	bool setWriterType(const int type);
	void setWriteHeader(const bool write){ this->output_header = write; }

	filter_type& getFilter(void){ return this->filter; }
	bool OpenWriter(void);
	bool OpenWriter(const std::string output_file);

private:
	bool __Open(const std::string input);
	bool ParseHeader(void);
	bool ParseHeaderExtend(void);
	bool __ParseRegion(const std::string& region, interval_type& interval);
	bool __ParseRegionIndexed(const std::string& region, interval_type& interval);
	bool __ParseRegionIndexedBlocks(void);
	bool __viewOnly(void);
	bool __viewFilter(void);
	bool __viewRegion(void);
	bool __viewRegionIndexed(void);
	bool __checkRegionIndex(const entry_type& entry);
	bool __checkRegionNoIndex(const entry_type& entry);
	bool __concat(const std::vector<std::string>& files, const std::string& output);

	bool addRegionsIndexed(std::vector<std::string>& positions);
	bool addRegionsUnindexed(std::vector<std::string>& positions);

public:
	U64 filesize;	// input file size
	U64 iterator_position_block;
	U64 iterator_position_variant;
	U64 size;
	bool hasIndex;
	std::ifstream stream; // reader stream
	header_type header; // header
	bool output_header;
	buffer_type compressed_buffer; // internal buffer
	buffer_type data_buffer; // internal buffer
	tgzf_controller_type tgzf_controller; // TGZF controller
	filter_type filter;	// filter parameters
	WRITER_TYPE writer_output_type;
	std::string literals; // header literals
	writer_type* writer; // writer interface
	contig_type* contigs;
	hash_table* contig_htable; // map input string to internal contigID
	tree_type** interval_tree;
	std::vector<interval_type>* interval_tree_entries;
	std::vector<totempole_sorted_entry_type>* interval_totempole_enties;
	toi_reader_type toi_reader;
};

}
} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTREADER_H_ */
