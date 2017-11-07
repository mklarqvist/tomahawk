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
	TomahawkImportEncoderStreamContainer() :
		stream_data_type(0),
		n_entries(0),
		n_additions(0)
	{}

	TomahawkImportEncoderStreamContainer(const U32 start_size) :
		stream_data_type(0),
		n_entries(0),
		n_additions(0),
		buffer(start_size)
	{}

	~TomahawkImportEncoderStreamContainer(){ this->buffer.deleteAll(); }

	inline void setType(const BYTE value){ this->stream_data_type = value; }

	inline void operator++(void){ ++this->n_additions; }

	inline void operator+=(const SBYTE& value){
		if(this->stream_data_type != 0) exit(1);
		this->buffer += value;
		++this->n_entries;
	}

	inline void operator+=(const BYTE& value){
		if(this->stream_data_type != 1) exit(1);
		this->buffer += value;
		++this->n_entries;
	}

	inline void operator+=(const S16& value){
		if(this->stream_data_type != 2) exit(1);
		this->buffer += value;
		++this->n_entries;
	}

	inline void operator+=(const U16& value){
		if(this->stream_data_type != 3) exit(1);
		this->buffer += value;
		++this->n_entries;
	}

	inline void operator+=(const S32& value){
		if(this->stream_data_type != 4) exit(1);
		this->buffer += value;
		++this->n_entries;
	}

	inline void operator+=(const U32& value){
		if(this->stream_data_type != 5) exit(1);
		this->buffer += value;
		++this->n_entries;
	}

	inline void operator+=(const U64& value){
		if(this->stream_data_type != 6) exit(1);
		this->buffer += value;
		++this->n_entries;
	}

	inline void operator+=(const float& value){
		if(this->stream_data_type != 7) exit(1);
		this->buffer += value;
		++this->n_entries;
	}

	void reset(void){
		this->stream_data_type = 0;
		this->n_entries = 0;
		this->n_additions = 0;
		this->buffer.reset();
	}

	inline void resize(const U32 size){ this->buffer.resize(size); }

public:
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
	typedef std::vector< std::vector<U32> > pattern_vector;

public:
	TomahawkImportWriter(const filter_type& filter);
	~TomahawkImportWriter();

	bool Open(const std::string output);
	void DetermineFlushLimit(void);
	bool OpenExtend(const std::string output);
	void WriteHeaders(void);
	void WriteFinal(void);
	void setHeader(vcf_header_type& header);
	bool add(const VCF::VCFLine& line);
	bool add(bcf_entry_type& entry, const U32* const ppa);
	bool parseBCF(meta_base_type& meta, bcf_entry_type& entry);

	void reset(void){
		// Buffers
		this->buffer_encode_rle.reset();
		this->buffer_encode_simple.reset();
		this->buffer_meta.reset();
		this->buffer_metaComplex.reset();

		// Hashes
		this->filter_hash_pattern_counter = 0;
		this->info_hash_pattern_counter = 0;
		this->format_hash_pattern_counter = 0;
		this->info_hash_value_counter = 0;
		this->format_hash_value_counter = 0;

		// Reset hash patterns
		this->filter_hash_pattern.clear();
		this->info_hash_pattern.clear();
		this->info_hash_streams.clear();
		this->format_hash_pattern.clear();
		this->format_hash_streams.clear();

		// Pattern data
		this->info_patterns.clear();
		this->format_patterns.clear();
		this->filter_patterns.clear();

		// Reset
		for(U32 i = 0; i < 100; ++i)
			this->containers[i].reset();
	}

	inline void TotempoleSwitch(const U32 contig, const U32 minPos){
		this->totempole_entry.reset();
		this->totempole_entry.contigID = contig;
		this->totempole_entry.minPosition = minPos;
	}

	// flush and write
	bool flush(const U32* const ppa);

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

	inline totempole_entry_type& getTotempoleEntry(void){ return(this->totempole_entry); }

public:
	// Stream information
	std::ofstream streamTomahawk;   // stream for data
	std::ofstream streamTotempole;  // stream for index
	std::string filename;
	std::string basePath;
	std::string baseName;

	// Basic
	U32 flush_limit;
	U32 n_variants_limit;
	U64 n_blocksWritten;            // number of blocks written
	U64 n_variants_written;         // number of variants written
	U64 n_variants_complex_written; // number of complex variants written
	U32 largest_uncompressed_block; // size of largest block observed
	const filter_type& filter;      // filters

	// Totempole entry
	totempole_entry_type totempole_entry;

	// TGZF controller
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
	// Hashes
	hash_table_vector filter_hash_pattern;
	hash_table_vector info_hash_pattern;
	hash_table info_hash_streams;
	hash_table_vector format_hash_pattern;
	hash_table format_hash_streams;
	// Actual patterns
	pattern_vector format_patterns;
	pattern_vector info_patterns;
	pattern_vector filter_patterns;

	// Encoder
	encoder_type* encoder;

	// VCF header reference
	vcf_header_type* vcf_header;

	// Stream containers
	stream_container* containers;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKIMPORTWRITER_H_ */
