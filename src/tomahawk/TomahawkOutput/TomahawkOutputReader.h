#ifndef TOMAHAWKOUTPUTREADER_H_
#define TOMAHAWKOUTPUTREADER_H_

#include <vector>
#include <fstream>
#include <string>
#include <stddef.h>
#include <regex>

#include "../../support/TypeDefinitions.h"
#include "../../io/BasicBuffer.h"
#include "../../io/compression/TGZFController.h"
#include "../../support/MagicConstants.h"
#include "../../algorithm/OpenHashTable.h"
#include "../../totempole/TotempoleMagic.h"
#include "../../third_party/intervalTree.h"
#include "../../totempole/TotempoleOutputReader.h"
#include "TomahawkOutputEntry.h"
#include "TomahawkOutputFilterController.h"
#include "TomahawkOutputWriter.h"
#include "../../totempole/TotempoleOutputIndexWriter.h"

namespace Tomahawk {
namespace IO {

class TomahawkOutputReader {
	typedef TomahawkOutputEntry entry_type;
	typedef TomahawkOutputFilterController filter_type;
	typedef Tomahawk::IO::TomahawkOutputWriterInterface writer_type;
	typedef TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> header_type;
	typedef Totempole::TotempoleContigBase contig_type;
	typedef TGZFHeader tgzf_type;
	typedef Hash::HashTable<std::string, U32> hash_table;
	typedef IO::TGZFController tgzf_controller_type;
	typedef Tomahawk::Algorithm::ContigInterval interval_type;
	typedef Tomahawk::Algorithm::IntervalTree<interval_type, U32> tree_type;
	typedef Totempole::TotempoleOutputSortedEntry totempole_sorted_entry_type;
	typedef IO::TomahawkOutputWriterIndex twoi_writer_type;
	typedef Tomahawk::IO::TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> toi_header_type;


public:
	typedef Totempole::TotempoleOutputReader toi_reader_type;
	enum WRITER_TYPE {binary, natural};

public:
	TomahawkOutputReader();
	~TomahawkOutputReader();

	const entry_type* operator[](const U32 p) const{ return(reinterpret_cast<const entry_type*>(&this->output_buffer.data[sizeof(entry_type)*p])); }

	// Streaming functions
	bool getBlock(const U32 blockID);
	bool getBlock(std::vector< std::pair<U32, U32> >& pairs);
	bool AddRegions(std::vector<std::string>& positions);
	bool Open(const std::string input);
	bool OpenExtend(const std::string input);
	bool nextBlock(const bool clear = true);
	bool nextVariant(const entry_type*& entry);
	bool nextVariantLimited(const entry_type*& entry);
	bool nextBlockUntil(const U32 limit);
	bool nextBlockUntil(const U32 limit, const U64 virtual_offset);
	inline void addLiteral(const std::string& string){ this->literals += string; }

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
	bool __checkRegionIndex(const entry_type* const entry);
	bool __checkRegionNoIndex(const entry_type* const entry);
	bool __concat(const std::vector<std::string>& files, const std::string& output);

	bool AddRegionsIndexed(std::vector<std::string>& positions);
	bool AddRegionsUnindexed(std::vector<std::string>& positions);

public:
	U64 filesize;	// input file size
	U64 position;
	U64 size;
	bool hasIndex;
	std::ifstream stream; // reader stream
	header_type header; // header
	bool output_header;
	IO::BasicBuffer buffer; // internal buffer
	IO::BasicBuffer output_buffer; // internal buffer
	tgzf_controller_type gzip_controller; // TGZF controller
	filter_type filter;	// filter parameters
	WRITER_TYPE writer_output_type;
	std::string literals; // header literals
	writer_type* writer; // writer interface
	contig_type* contigs;
	hash_table* contig_htable; // map input string to internal contigID
	tree_type** interval_tree;
	std::vector<interval_type>* interval_tree_entries;
	std::vector<totempole_sorted_entry_type>* interval_totempole_entries;
	toi_reader_type toi_reader;
};

}
} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTREADER_H_ */
