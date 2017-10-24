#ifndef TOMAHAWKIMPORTWRITER_H_
#define TOMAHAWKIMPORTWRITER_H_

#include <fstream>

#include "../algorithm/compression/TomahawkImportEncoder.h"
#include "../support/TypeDefinitions.h"
#include "../io/BasicBuffer.h"
#include "../io/BasicWriters.h"
#include "../io/compression/TGZFController.h"
#include "../io/vcf/VCFHeaderConstants.h"
#include "../io/vcf/VCFLines.h"
#include "../io/vcf/VCFHeader.h"
#include "base/TomahawkEntryMeta.h"
#include "base/TomahawkEntrySupport.h"
#include "../totempole/TotempoleEntry.h"
#include "../totempole/TotempoleReader.h"
#include "../support/simd_definitions.h"
#include "TomahawkImporterFilters.h"

namespace Tomahawk {

class TomahawkImportWriter {
	typedef IO::BasicBuffer buffer_type;
	typedef TomahawkImporterFilters filter_type;
	typedef Support::TomahawkEntryMetaBase meta_base_type;
	typedef Algorithm::TomahawkImportEncoder encoder_type;
	typedef VCF::VCFHeader vcf_header_type;
	typedef Totempole::TotempoleEntry totempole_entry_type;
	typedef IO::TGZFController tgzf_controller_type;
	typedef BCF::BCFEntry bcf_entry_type;

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
	bool add(const bcf_entry_type& line);

	inline void reset(void){
		this->buffer_encode_rle.reset();
		this->buffer_encode_simple.reset();
		this->buffer_meta.reset();
		this->buffer_metaComplex.reset();
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
		if(this->totempole_entry.n_variantsRLE + this->totempole_entry.n_variantsComplex >= this->n_variants_limit || this->buffer_encode_rle.size() >= this->flush_limit){
			//std::cerr << "flushing: " << this->totempole_entry_.variants << '/' << this->n_variants_limit << '\t' << this->buffer_rle.size() << '/' << this->flush_limit << std::endl;
			return true;
		}

		return false;
	}

	inline const U64& blocksWritten(void) const{ return this->n_blocksWritten; }
	inline const U64& size(void) const{ return this->buffer_encode_rle.size(); }
	inline const U64& getVariantsWritten(void) const{ return this->n_variants_written; }

	void CheckOutputNames(const std::string& input);

	inline Totempole::TotempoleEntry& getTotempoleEntry(void){ return(this->totempole_entry); }

public:
	std::ofstream streamTomahawk;   // stream for data
	std::ofstream streamTotempole;  // stream for index
	std::string filename;
	std::string basePath;
	std::string baseName;

	U32 flush_limit;
	U32 n_variants_limit;
	U64 n_blocksWritten;            // number of blocks written
	U64 n_variants_written;         // number of variants written
	U64 n_variants_complex_written; // number of complex variants written
	U32 largest_uncompressed_block; // size of largest block in b
	const filter_type& filter;      // filters

	totempole_entry_type totempole_entry;
	tgzf_controller_type gzip_controller;

	buffer_type buffer_encode_rle; // run lengths
	buffer_type buffer_encode_simple; // simple encoding
	buffer_type buffer_meta;// meta data for run lengths (chromosome, position, ref/alt)
	buffer_type buffer_metaComplex; // complex meta data

	encoder_type* encoder;
	vcf_header_type* vcf_header;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKIMPORTWRITER_H_ */
