#ifndef TOMAHAWKIMPORTWRITER_H_
#define TOMAHAWKIMPORTWRITER_H_

#include <fstream>

#include "../support/TypeDefinitions.h"
#include "../io/BasicBuffer.h"
#include "../io/BasicWriters.h"
#include "../io/compression/TGZFController.h"
#include "../io/vcf/VCFHeaderConstants.h"
#include "../io/vcf/VCFLines.h"
#include "../io/vcf/VCFHeader.h"
#include "../algorithm/compression/TomahawkImportRLE.h"
#include "base/TomahawkEntryMeta.h"
#include "../totempole/TotempoleEntry.h"
#include "../totempole/TotempoleReader.h"
#include "../support/simd_definitions.h"
#include "TomahawkImporterFilters.h"

namespace Tomahawk {

class TomahawkImportWriter {
	typedef IO::BasicBuffer buffer_type;
	typedef TomahawkImporterFilters filter_type;
	typedef Support::TomahawkEntryMetaBase meta_base_type;

public:
	TomahawkImportWriter(const filter_type& filter);
	~TomahawkImportWriter();

	bool Open(const std::string output);
	void DetermineFlushLimit(void);
	bool OpenExtend(const std::string output);
	void WriteHeaders(void);
	void WriteFinal(void);
	void setHeader(VCF::VCFHeader& header);
	bool add(const VCF::VCFLine& line);
	bool add(const BCF::BCFEntry& line);

	inline void reset(void){
		this->buffer_rle_.reset();
		this->buffer_meta_.reset();
	}

	inline void TotempoleSwitch(const U32 contig, const U32 minPos){
		this->totempole_entry.reset();
		this->totempole_entry.contigID = contig;
		this->totempole_entry.minPosition = minPos;
	}

	// flush and write
	bool flush(void);

	inline bool checkSize() const{
		// if the current size is larger than our desired output block size, return TRUE to trigger a flush
		// or if the number of entries written to buffer exceeds our set limit
		if(this->totempole_entry.variants >= this->n_variants_limit || this->buffer_rle_.size() >= this->flush_limit){
			//std::cerr << "flushing: " << this->totempole_entry_.variants << '/' << this->n_variants_limit << '\t' << this->buffer_rle_.size() << '/' << this->flush_limit << std::endl;
			return true;
		}

		return false;
	}

	inline const U64& blocksWritten(void) const{ return this->blocksWritten_; }
	inline const U64& size(void) const{ return this->buffer_rle_.size(); }
	inline const U64& getVariantsWritten(void) const{ return this->variants_written; }

	void CheckOutputNames(const std::string& input);

	inline Totempole::TotempoleEntry& getTotempoleEntry(void){ return(this->totempole_entry); }

public:
	std::ofstream streamTomahawk;	// stream
	std::ofstream streamTotempole;	// stream
	U32 flush_limit;
	U32 n_variants_limit;
	U64 blocksWritten_;				// number of blocks written
	U64 variants_written;			// number of variants written
	U64 variants_complex_written;	// number of complex variants written
	U32 largest_uncompressed_block_;// size of largest block in b
	const filter_type& filter;		// filters

	Totempole::TotempoleEntry totempole_entry;
	IO::TGZFController gzip_controller_;
	Algorithm::TomahawkImportRLE* rleController_;
	IO::BasicBuffer buffer_rle_;	// run lengths
	IO::BasicBuffer buffer_meta_;	// meta data for run lengths (chromosome, position, ref/alt)

	VCF::VCFHeader* vcf_header_;

	std::string filename;
	std::string basePath;
	std::string baseName;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKIMPORTWRITER_H_ */
