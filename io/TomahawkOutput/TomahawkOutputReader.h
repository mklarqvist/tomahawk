#ifndef TOMAHAWKOUTPUTREADER_H_
#define TOMAHAWKOUTPUTREADER_H_

#include <vector>
#include <fstream>
#include <string>
#include <stddef.h>
#include <regex>

#include "../../TypeDefinitions.h"
#include "../../tomahawk/MagicConstants.h"
#include "TomahawkOutputEntry.h"
#include "../../io/PackedEntryReader.h"
#include "TomahawkOutputFilterController.h"
#include "../../io/BasicBuffer.h"
#include "../../io/GZController.h"
#include "../../io/TomahawkOutput/TomahawkOutputWriter.h"
#include "../../totempole/TotempoleMagic.h"
#include "../../io/TGZFHeader.h"
#include "../../third_party/intervalTree.h"


namespace Tomahawk {
namespace IO {
// Todo: TomahawkOutputIndexReader

class TomahawkOutputReader {
	typedef TomahawkOutputEntry entry_type;
	typedef TomahawkOutputFilterController filter_type;
	typedef PackedEntryReader<entry_type, sizeof(entry_type)> reader_type;
	typedef Tomahawk::IO::TomahawkOutputWriterInterface writer_type;
	typedef TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> header_type;
	typedef Totempole::TotempoleContigBase contig_type;
	typedef TGZFHeader tgzf_type;
	typedef Hash::HashTable<std::string, U32> hash_table;
	typedef IO::GZController tgzf_controller_type;
	typedef Tomahawk::Algorithm::ContigInterval interval_type;
	typedef Tomahawk::Algorithm::IntervalTree<interval_type, U32> tree_type;

	enum WRITER_TYPE {binary, natural};

public:
	TomahawkOutputReader();
	~TomahawkOutputReader();

	const entry_type* operator[](const U32 p) const{ return(reinterpret_cast<const entry_type*>(&this->output_buffer.data[sizeof(entry_type)*p])); }

	// Streaming functions
	bool getBlock(const U32 blockID);
	bool getBlock(std::vector< std::pair<U32, U32> >& pairs);
	bool getBlocks(void);
	bool AddRegions(std::vector<std::string>& positions);

	bool Open(const std::string input);
	bool nextBlock(void);
	bool nextVariant(const entry_type*& entry);

	// Other
	bool view(const std::string& filename);
	bool index(const std::string& filename);
	bool summary(const std::string& input);

	//
	bool setWriterType(const int type){
		if(type == 0)
			this->writer_output_type = WRITER_TYPE::binary;
		else if(type == 1)
			this->writer_output_type = WRITER_TYPE::natural;
		else {
			std::cerr << "unknown writer type: " << type << std::endl;
			return false;
		}
		return true;
	}
	void setWriteHeader(const bool write){ this->output_header = write; }

	// Read entire file into memory
	filter_type& getFilter(void){ return this->filter; }

private:
	bool ParseHeader(void);
	bool __ParseRegion(const std::string& region, interval_type& interval);
	bool __viewOnly(void);
	bool __viewFilter(void);
	bool __viewRegion(void);
	bool __checkRegion(const entry_type* const entry);
	void __openWriter(void);

public:
	U64 filesize;	// input file size
	U32 position;
	U32 size;

	std::ifstream stream; // reader stream
	reader_type reader; // reader
	header_type header; // header
	bool output_header;

	IO::BasicBuffer buffer; // internal buffer
	IO::BasicBuffer output_buffer; // internal buffer
	tgzf_controller_type gzip_controller; // TGZF controller
	filter_type filter;	// filter parameters
	WRITER_TYPE writer_output_type;
	writer_type* writer; // writer interface
	// Todo: PackedEntryIterator taking as input char* and length or IO::BasicBuffer
	contig_type* contigs;
	hash_table* contig_htable; // map input string to internal contigID
	tree_type** interval_tree;
	std::vector<interval_type>* interval_tree_entries;
};

}
} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTREADER_H_ */
