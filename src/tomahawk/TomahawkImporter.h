#ifndef VCFPARSER_H_
#define VCFPARSER_H_

#include "../io/reader.h"
#include "../io/vcf/VCFHeader.h"
#include "../io/bcf/BCFReader.h"
#include "TomahawkImporterFilters.h"
#include "TomahawkImportWriter.h"

namespace Tomahawk {

class TomahawkImporter {
	typedef TomahawkImporter self_type;
	typedef reader reader_type;
	typedef VCF::VCFHeader header_type;
	typedef TomahawkImportWriter writer_type;
	typedef VCF::VCFLine line_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Algorithm::TomahawkImportRLE rle_controller_type;
	typedef Totempole::TotempoleEntry totempole_entry_type;
	typedef BCF::BCFReader bcf_reader_type;
	typedef BCF::BCFEntry bcf_entry_type;
	typedef TomahawkImporterFilters filter_type;

	/*
	 This supportive structure keeps track of the current and
	 previous contig identifiers and the previous obseved position.
	 This information is necessary to guarantee the sort-order of
	 the output Tomahawk file required for indexing.
	 The flag previous_included is triggered whenever an entry is
	 not filtered out. It is used when two or more entries share the
	 same position. In this case, if the preceding line was included
	 then ignore the current one. Otherwise, the preceding line was
	 filtered out and the include the current one.
	 Note that contigID is a pointer as this is required by our
	 hash-table implementation as a return value
	 */
	struct __InternalHelper{
		__InternalHelper():
			contigID(nullptr),
			prevcontigID(-1),
			previous_position(-1),
			previous_included(false)
		{}
		S32* contigID;			// current contigID
		S32 prevcontigID;		// previous contigID
		S32 previous_position;	// current position
		bool previous_included;
	} sort_order_helper;

public:
	TomahawkImporter(std::string inputFile, std::string outputPrefix);
	~TomahawkImporter();
	bool Build();
	bool Extend(std::string extendFile);
	filter_type& getFilters(void){ return(this->filters); }

private:
	bool BuildVCF();  // import a VCF file
	bool BuildBCF();  // import a BCF file
	bool ExtendVCF(); // extend a Twk file with a VCF file
	bool ExtendBCF(); // extend a Twk file with a BCF file

	bool parseVCFLine(line_type& line); // Import a VCF line
	bool parseBCFLine(bcf_entry_type& line); // Import a BCF line
	// Check if the current meta and RLE buffers exceeds
	// the disk flush limit
	bool checkSize(void) const{ return(this->meta_buffer.size() + this->rle_buffer.size() >= this->block_flush_limit); }

private:
	U32 block_flush_limit;    // limit in bytes when to flush to disk
	std::string inputFile;    // input file name
	std::string outputPrefix; // output file prefix
	reader_type reader_;      // reader
	writer_type writer_;      // writer
	buffer_type meta_buffer;  // meta buffer
	buffer_type rle_buffer;   // RLE buffer
	totempole_entry_type totempole_entry;  // totempole entry for indexing
	filter_type filters;
	header_type* header_;     // header
	rle_controller_type* rle_controller;   // RLE packer
};


} /* namespace Tomahawk */

#endif /* VCFPARSER_H_ */
