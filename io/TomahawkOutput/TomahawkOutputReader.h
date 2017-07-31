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
#include "../../io/BasicWriters.h"
#include "../../io/totempole/TotempoleMagic.h"
#include "../../io/TGZFHeader.h"
#include "../../third_party/intervalTree.h"


namespace Tomahawk {
namespace IO {
// Todo: TomahawkOutputIndexReader

class TomahawkOutputReader {
	typedef TomahawkOutputEntry entry_type;
	typedef TomahawkOutputFilterController filter_type;
	typedef PackedEntryReader<entry_type, sizeof(entry_type)> reader_type;
	typedef IO::GenericWriterInterace writer_type;
	typedef TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> header_type;
	typedef Totempole::TotempoleContigBase contig_type;
	typedef TGZFHeader tgzf_type;
	typedef Hash::HashTable<std::string, U32> hash_table;
	typedef IO::GZController tgzf_controller_type;
	typedef Tomahawk::Algorithm::ContigInterval interval_type;
	typedef Tomahawk::Algorithm::IntervalTree<interval_type, U32> tree_type;

public:
	TomahawkOutputReader();
	~TomahawkOutputReader(){
		delete [] this->contigs;
		delete contig_htable;
		if(interval_tree != nullptr){
			for(U32 i = 0; i < this->header.n_contig; ++i)
				delete this->interval_tree[i];
		}
		delete [] interval_tree_entries;

		delete interval_tree;
		delete writer;

		this->buffer.deleteAll();
		this->output_buffer.deleteAll();
	}

	const TomahawkOutputEntry* operator[](const U32 p) const{ return(reinterpret_cast<TomahawkOutputEntry*>(&this->output_buffer.data[sizeof(TomahawkOutputEntry)*p])); }

	// Streaming functions
	bool getBlock(const U32 blockID);
	bool getBlock(std::vector< std::pair<U32, U32> >& pairs);
	bool getBlocks(void);
	bool AddRegions(std::vector<std::string>& positions);

	bool Open(const std::string input);
	bool ParseHeader(void);
	bool nextBlock(void);
	bool nextVariant(const TomahawkOutputEntry*& entry);

	// Other
	bool view(const std::string& filename);
	bool index(const std::string& filename);
	bool summary(const std::string& input);

	// Read entire file into memory
	filter_type& getFilter(void){ return this->filter; }

private:
	bool __ParseRegion(const std::string& region, interval_type& interval);
	bool __viewOnly(void);
	bool __viewFilter(void);
	bool __viewRegion(void);

public:
	//U64 samples; 	// has to match header
	//float version;	// has to match header
	U64 filesize;	// input file size

	U32 position;
	U32 size;

	std::ifstream stream; // reader stream
	header_type header; // header

	IO::BasicBuffer buffer; // internal buffer
	IO::BasicBuffer output_buffer; // internal buffer
	tgzf_controller_type gzip_controller; // TGZF controller
	filter_type filter;	// filter parameters
	writer_type* writer; // writer interface
	// Todo: PackedEntryIterator taking as input char* and length or IO::BasicBuffer
	contig_type* contigs;
	hash_table* contig_htable; // map input string to internal contigID
	tree_type** interval_tree;
	std::vector<interval_type>* interval_tree_entries;

	//temp
	reader_type reader;
};

}
} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTREADER_H_ */
