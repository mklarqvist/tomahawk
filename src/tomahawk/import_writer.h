#ifndef TOMAHAWKIMPORTWRITER_H_
#define TOMAHAWKIMPORTWRITER_H_

#include <fstream>

#include "../algorithm/compression/genotype_encoder.h"
#include "../support/type_definitions.h"
#include "../io/BasicBuffer.h"
#include "../io/BasicWriters.h"
#include "../io/compression/TGZFController.h"
#include "../io/vcf/VCFHeaderConstants.h"
#include "../io/vcf/VCFLines.h"
#include "../io/vcf/VCFHeader.h"
#include "../index/index_entry.h"
#include "../totempole/TotempoleReader.h"
#include "../support/simd_definitions.h"
#include "import_filters.h"
#include "meta_entry.h"

namespace Tomahawk {

class ImportWriter {
private:
	typedef IO::BasicBuffer       buffer_type;
	typedef ImporterFilters       filter_type;
	typedef Totempole::IndexEntry index_entry_type;

public:
	ImportWriter(const filter_type& filter);
	~ImportWriter();

	bool Open(const std::string output);
	void DetermineFlushLimit(void);
	bool OpenExtend(const std::string output);
	int WriteHeaders(void);
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
		this->totempole_entry.min_position = minPos;
	}

	// flush and write
	bool flush(void);

	inline bool checkSize() const{
		// if the current size is larger than our desired output block size, return TRUE to trigger a flush
		// or if the number of entries written to buffer exceeds our set limit
		if(this->totempole_entry.n_variants >= this->n_variants_limit || this->buffer_rle_.size() >= this->flush_limit){
			//std::cerr << "flushing: " << this->totempole_entry_.variants << '/' << this->n_variants_limit << '\t' << this->buffer_rle_.size() << '/' << this->flush_limit << std::endl;
			return true;
		}

		return false;
	}

	inline const U32& blocksWritten(void) const{ return this->blocksWritten_; }
	inline const U64& size(void) const{ return this->buffer_rle_.size(); }

	void CheckOutputNames(const std::string& input);


	inline U32 GetVariantsWritten(void) const{ return this->variants_written_; }
	inline index_entry_type& getTotempoleEntry(void){ return(this->totempole_entry); }

public:
	std::ofstream stream;	// stream

	U32 flush_limit;
	U32 n_variants_limit;
	U32 blocksWritten_;				// number of blocks written
	U32 variants_written_;			// number of variants written
	U32 largest_uncompressed_block_;// size of largest block in b
	const filter_type& filter;		// filters

	index_entry_type totempole_entry;
	IO::TGZFController gzip_controller_;
	Algorithm::GenotypeEncoder* rleController_;
	buffer_type buffer_rle_;	// run lengths
	buffer_type buffer_meta_;	// meta data for run lengths (chromosome, position, ref/alt)

	VCF::VCFHeader* vcf_header_;

	std::string filename;
	std::string basePath;
	std::string baseName;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKIMPORTWRITER_H_ */
