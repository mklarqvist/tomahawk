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
#include "../algorithm/OpenHashTable.h"

namespace Tomahawk {

// Stream container for importing
class TomahawkImportEncoderStreamContainer{
	typedef TomahawkImportEncoderStreamContainer self_type;
	typedef IO::BasicBuffer buffer_type;

public:
	TomahawkImportEncoderStreamContainer(const U32 start_size) :
		stream_data_type(0),
		n_entries(0),
		n_additions(0),
		buffer(start_size)
	{}
	~TomahawkImportEncoderStreamContainer(){ this->buffer.deleteAll(); }

	template <class T>
	inline void operator+=(const T& value){ this->buffer += value; }

	void reset(void){
		this->stream_data_type = 0;
		this->n_entries = 0;
		this->n_additions = 0;
		this->buffer.reset();
	}

private:
	BYTE stream_data_type;
	U32 n_entries;   // number of entries
	U32 n_additions; // number of times added
	buffer_type buffer;
};

class TomahawkImportWriter {
	typedef IO::BasicBuffer buffer_type;
	typedef TomahawkImporterFilters filter_type;
	typedef Support::TomahawkEntryMetaBase meta_base_type;
	typedef Algorithm::TomahawkImportEncoder encoder_type;
	typedef VCF::VCFHeader vcf_header_type;
	typedef Totempole::TotempoleEntry totempole_entry_type;
	typedef IO::TGZFController tgzf_controller_type;
	typedef BCF::BCFEntry bcf_entry_type;
	typedef Hash::HashTable<U64, U32> hash_table_vector;
	typedef Hash::HashTable<U32, U32> hash_table;
	typedef TomahawkImportEncoderStreamContainer stream_container;

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
	bool add(bcf_entry_type& entry);
	bool parseBCF(bcf_entry_type& entry);

	inline void reset(void){
		this->buffer_encode_rle.reset();
		this->buffer_encode_simple.reset();
		this->buffer_meta.reset();
		this->buffer_metaComplex.reset();

		this->filter_hash_pattern_counter = 0;
		this->info_hash_pattern_counter = 0;
		this->format_hash_pattern_counter = 0;
		this->info_hash_value_counter = 0;
		this->format_hash_value_counter = 0;
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
		if(this->totempole_entry.n_variants >= this->n_variants_limit || this->buffer_encode_rle.size() >= this->flush_limit)
			return true;

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
	U32 largest_uncompressed_block; // size of largest block observed
	const filter_type& filter;      // filters

	totempole_entry_type totempole_entry;
	tgzf_controller_type gzip_controller;

	// Basic output buffers
	buffer_type buffer_encode_rle; // run lengths
	buffer_type buffer_encode_simple; // simple encoding
	buffer_type buffer_meta;// meta data for run lengths (chromosome, position, ref/alt)
	buffer_type buffer_metaComplex; // complex meta data

	// Meta data
	U32 filter_hash_pattern_counter;
	U32 info_hash_pattern_counter;
	U32 format_hash_pattern_counter;
	U32 info_hash_value_counter;
	U32 format_hash_value_counter;
	hash_table_vector filter_hash_pattern;
	hash_table_vector info_hash_pattern;
	hash_table info_hash_streams;
	hash_table_vector format_hash_pattern;
	hash_table format_hash_streams;

	encoder_type* encoder;
	vcf_header_type* vcf_header;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKIMPORTWRITER_H_ */
