#ifndef TOMAHAWK_TOMAHAWK_IMPORTER_
#define TOMAHAWK_TOMAHAWK_IMPORTER_

#include "io/bcf/bcf_reader.h"
#include "io/vcf/vcf_header.h"
#include "io/reader.h"
#include "import_filters.h"
#include "import_writer.h"
#include "index/index.h"

namespace tomahawk {

struct TomahawkImporterStats{
	TomahawkImporterStats() :
		n_filtered_eov(0),
		n_filtered_missingness(0),
		n_filtered_not_biallelic(0)
	{

	}

	~TomahawkImporterStats() = default;

	U64 n_filtered_eov;
	U64 n_filtered_missingness;
	U64 n_filtered_not_biallelic;
};

/**<
 * This class handles importing `bcf`/`vcf` into the `twk` file format.
 */
class TomahawkImporter {
	typedef TomahawkImporter           self_type;
	typedef reader                     reader_type;
	typedef ImportWriter               writer_type;
	typedef ImporterFilters            filter_type;
	typedef io::BasicBuffer            buffer_type;
	typedef algorithm::GenotypeEncoder rle_controller_type;
	typedef totempole::IndexEntry      totempole_entry_type;
	typedef vcf::VCFHeader             vcf_header_type;
	typedef vcf::VCFLine               vcf_entry_type;
	typedef bcf::BCFReader             bcf_reader_type;
	typedef bcf::BCFEntry              bcf_entry_type;
	typedef Index                      index_type;
	typedef totempole::Footer          footer_type;

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
	TomahawkImporter(const std::string inputFile, const std::string outputPrefix);
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
	bool parseVCFLine(vcf_entry_type& line); // Import a VCF line
	bool parseBCFLine(bcf_entry_type& line); // Import a BCF line

	// Check if the current meta and RLE buffers exceeds
	// the disk flush limit
	inline bool checkSize(void) const{ return(this->meta_buffer.size() + this->rle_buffer.size() >= this->block_flush_limit); }

private:
	U32                  block_flush_limit;// limit in bytes when to flush to disk
	std::string          input_file;        // input file name
	std::string          output_prefix;     // output file prefix
	reader_type          reader_;          // reader
	writer_type          writer_;          // writer
	buffer_type          meta_buffer;      // meta buffer
	buffer_type          rle_buffer;       // RLE buffer
	totempole_entry_type totempole_entry;  // current (active) index entry
	filter_type          filters;          // filters
	index_type           index;
	footer_type          footer_;
	vcf_header_type*     vcf_header_;      // vcf header
	//rle_controller_type* rle_controller;   // RLE packer algorithms
};


} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWK_IMPORTER_ */
