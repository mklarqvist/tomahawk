#ifndef VCFPARSER_H_
#define VCFPARSER_H_

#include "VCFHeader.h"
#include "../../tomahawk/TomahawkImportWriter.h"

namespace Tomahawk {
namespace VCF{

class VCFParser {
	typedef VCFParser self_type;
	typedef reader reader_type;
	typedef VCFHeader header_type;
	typedef TomahawkImportWriter writer_type;
	typedef VCFLine line_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Algorithm::TomahawkImportRLE rle_controller_type;
	typedef TotempoleEntry totempole_entry_type;

	struct __InternalHelper{
		__InternalHelper(): contigID(nullptr), prevcontigID(nullptr), previous_position(0){}
		U32* contigID;			// current contigID
		U32* prevcontigID;		// previous contigID
		U32 previous_position;	// current position
	} sort_order_helper;

public:
	VCFParser(std::string inputFile, std::string outputPrefix);
	~VCFParser();
	bool Build();
	bool Extend();

private:
	bool parseLine(line_type& line);
	bool checkSize(void) const{ return(this->meta_buffer.size() + this->rle_buffer.size() >= this->block_flush_limit); }

private:
	U32 block_flush_limit;
	std::string inputFile;
	std::string outputPrefix;
	reader_type reader_;
	header_type header_;
	writer_type writer_;

	buffer_type meta_buffer;
	buffer_type rle_buffer;
	totempole_entry_type totempole_entry;
	rle_controller_type* rle_controller;
};


} /* namespace VCF */
} /* namespace Tomahawk */

#endif /* VCFPARSER_H_ */
