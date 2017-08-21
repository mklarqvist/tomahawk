#ifndef VCFPARSER_H_
#define VCFPARSER_H_

#include "../io/vcf/VCFHeader.h"
#include "../io/bcf/BCFReader.h"
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
	typedef TotempoleEntry totempole_entry_type;
	typedef BCF::BCFReader bcf_reader_type;
	typedef BCF::BCFEntry bcf_entry_type;

	struct __InternalHelper{
		__InternalHelper(): contigID(nullptr), prevcontigID(-1), previous_position(-1), previous_included(false){}
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

private:
	bool BuildVCF();
	bool BuildBCF();
	bool ExtendVCF();
	bool ExtendBCF();

	bool parseVCFLine(line_type& line);
	bool parseBCFLine(bcf_entry_type& line);
	bool checkSize(void) const{ return(this->meta_buffer.size() + this->rle_buffer.size() >= this->block_flush_limit); }

private:
	U32 block_flush_limit;
	std::string inputFile;
	std::string outputPrefix;
	reader_type reader_;
	header_type* header_;
	writer_type writer_;

	buffer_type meta_buffer;
	buffer_type rle_buffer;
	totempole_entry_type totempole_entry;
	rle_controller_type* rle_controller;
};


} /* namespace Tomahawk */

#endif /* VCFPARSER_H_ */
