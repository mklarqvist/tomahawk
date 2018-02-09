#ifndef VCFPARSER_H_
#define VCFPARSER_H_

#include "../io/reader.h"
#include "../io/vcf/VCFHeader.h"
#include "../io/bcf/BCFReader.h"
#include "import_filters.h"
#include "import_writer.h"

namespace Tomahawk {

/**<
 * This class handles importing `bcf`/`vcf` into `twk` format.
 */
class TomahawkImporter {
	typedef TomahawkImporter           self_type;
	typedef reader                     reader_type;
	typedef VCF::VCFHeader             header_type;
	typedef VCF::VCFLine               line_type;
	typedef ImportWriter               writer_type;
	typedef IO::BasicBuffer            buffer_type;
	typedef Algorithm::GenotypeEncoder rle_controller_type;
	typedef Totempole::IndexEntry  totempole_entry_type;
	typedef BCF::BCFReader             bcf_reader_type;
	typedef BCF::BCFEntry              bcf_entry_type;
	typedef ImporterFilters            filter_type;

	/**<
	 * This supportive structure keeps track of the current and
	 * previous contig identifiers and the previous obseved position.
	 * This information is necessary to guarantee the sort-order of
	 * the output Tomahawk file required for indexing.
	 * The flag previous_included is triggered whenever an entry is
	 * not filtered out. It is used when two or more entries share the
	 * same position. In this case, if the preceding line was included
	 * then ignore the current one. Otherwise, the preceding line was
	 * filtered out and the include the current one.
	 * Note that contigID is a pointer as this is required by our
	 * hash-table implementation as a return value
	 */
	struct __InternalHelper {
		__InternalHelper():
			contigID(nullptr),
			prevcontigID(-1),
			previous_position(-1),
			previous_included(false)
		{}
		S32* contigID;          // current contigID
		S32 prevcontigID;       // previous contigID
		S32 previous_position;  // current position
		bool previous_included;
	} sort_order_helper;

public:
	TomahawkImporter(std::string inputFile, std::string outputPrefix);
	~TomahawkImporter();

	/**<
	 * Primary import function for data. The function internally checks
	 * the target input file type (`bcf`/`vcf`).
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool Build();

	/**<
	 * Extends an existing `twk` file with a target input file.
	 * Warning: this function does NOT check if the file headers
	 * are in order! If they are not in order the `twk` file will
	 * be corrupted
	 * @param extendFile Target `twk` file to extend
	 * @return           Returns TRUE upon success or FALSE otherwise
	 */
	bool Extend(std::string extendFile);
	filter_type& getFilters(void){ return(this->filters); }

private:
	// Basic import funtionality
	bool BuildVCF();  // import a VCF file
	bool BuildBCF();  // import a BCF file

	// Extend existing `twk` file with data from a `vcf`/`bcf` file
	bool ExtendVCF(); // extend a Twk file with a VCF file
	bool ExtendBCF(); // extend a Twk file with a BCF file

	// Parse a `bcf`/`vcf` line
	bool parseVCFLine(line_type& line);      // Import a VCF line
	bool parseBCFLine(bcf_entry_type& line); // Import a BCF line

	// Check if the current meta and RLE buffers exceeds
	// the disk flush limit
	inline bool checkSize(void) const{ return(this->meta_buffer.size() + this->rle_buffer.size() >= this->block_flush_limit); }

private:
	U32 block_flush_limit;    // limit in bytes when to flush to disk
	std::string inputFile;    // input file name
	std::string outputPrefix; // output file prefix
	reader_type reader_;      // reader
	writer_type writer_;      // writer
	buffer_type meta_buffer;  // meta buffer
	buffer_type rle_buffer;   // RLE buffer
	totempole_entry_type totempole_entry;  // totempole entry for indexing
	filter_type filters;      // filters
	header_type* header_;     // header
	rle_controller_type* rle_controller; // RLE packer
};


} /* namespace Tomahawk */

#endif /* VCFPARSER_H_ */
