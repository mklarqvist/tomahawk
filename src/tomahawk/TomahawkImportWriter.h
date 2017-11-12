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
#include "../algorithm/permutation/RadixSortGT.h"

namespace Tomahawk {

inline bool bytePreprocessor(const U32* const data, const size_t& size, char* destination){
	if(size == 0) return false;

	BYTE* targetData = reinterpret_cast<BYTE*>(destination);
	BYTE* s1 = &targetData[size*3];
	BYTE* s2 = &targetData[size*2];
	BYTE* s3 = &targetData[size*1];
	BYTE* s4 = &targetData[size*0];

	for(U32 i = 0; i < size; ++i){
		const BYTE* const p = reinterpret_cast<const BYTE* const>(&data[i]);
		s1[i] = (p[0] & 255);
		s2[i] = (p[1] & 255);
		s3[i] = (p[2] & 255);
		s4[i] = (p[3] & 255);
		const U32 x = s4[i] << 24 | s3[i] << 16 | s2[i] << 8 | s1[i];
		assert(x == data[i]);
	}

	return true;
}

inline bool bytePreprocessorRevert(const char* const data, const size_t& size, char* destination){
	if(size == 0) return false;

	const BYTE* d = reinterpret_cast<const BYTE*>(data);
	U32* dest = reinterpret_cast<U32*>(destination);
	const BYTE* const s1 = &d[size*3];
	const BYTE* const s2 = &d[size*2];
	const BYTE* const s3 = &d[size*1];
	const BYTE* const s4 = &d[size*0];

	for(U32 i = 0; i < size; ++i){
		const U32 x = s4[i] << 24 | s3[i] << 16 | s2[i] << 8 | s1[i];
		dest[i] = x;
	}

	return true;
}

inline bool bytePreprocessBits(const char* const data, const size_t& size, char* destination){
	if(size == 0) return false;

	BYTE* dest = reinterpret_cast<BYTE*>(destination);
	const BYTE* const d = reinterpret_cast<const BYTE* const>(data);
	BYTE* target[8];
	const U32 s = size/8;
	for(U32 i = 0; i < 8; ++i)
		target[7-i] = &dest[s*i];

	U32 k = 0; U32 p = 0;
	for(U32 i = 0; i < size; ++i, ++k){
		for(U32 j = 0; j < 8; ++j)
			target[j][p] |= ((d[i] & (1 << j)) >> j) << k;

		if(i % 7 == 0){ k = 0; ++p; }
	}

	return true;
}

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
	typedef std::vector<U32> id_vector;
	typedef std::vector< id_vector > pattern_vector;


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

		// Value data
		this->info_values.clear();
		this->format_values.clear();
		this->filter_values.clear();

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
	bool flush(Algorithm::RadixSortGT& permuter);
	//bool flush(const U32* const ppa);

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

private:
	bool constructBitVector(const hash_table_vector& htable, const U32 byte_width);

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
	buffer_type buffer_ppa; // ppa buffer

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
	// Values
	id_vector format_values;
	id_vector info_values;
	id_vector filter_values;

	// Encoder
	encoder_type* encoder;

	// VCF header reference
	vcf_header_type* vcf_header;

	// Stream containers
	stream_container* containers;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKIMPORTWRITER_H_ */
