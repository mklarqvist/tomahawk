#include "TomahawkImportWriter.h"

#include <strings.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "../io/bcf/BCFEntry.h"
#include "../support/helpers.h"
#include "../support/MagicConstants.h"
#include "../totempole/TotempoleContig.h"
#include "../totempole/TotempoleHeader.h"
#include "base/TomahawkGTEntries.h"

#include "TomahawkBlockIterator.h"

#include "../third_party/FiniteStateEntropy/lib/fse.h"

namespace Tomahawk {

TomahawkImportWriter::TomahawkImportWriter(const filter_type& filter) :
	flush_limit(1000000),
	n_variants_limit(1024),
	n_blocksWritten(0),
	n_variants_written(0),
	n_variants_complex_written(0),
	largest_uncompressed_block(0),
	filter(filter),
	buffer_encode_rle(flush_limit*2),
	buffer_encode_simple(flush_limit*2),
	buffer_meta(flush_limit*10), // meta joins all other buffers
	buffer_metaComplex(flush_limit*2),
	filter_hash_pattern_counter(0),
	info_hash_pattern_counter(0),
	format_hash_pattern_counter(0),
	info_hash_value_counter(0),
	format_hash_value_counter(0),
	filter_hash_value_counter(0),
	filter_hash_pattern(5012),
	info_hash_pattern(5012),
	info_hash_streams(5012),
	format_hash_pattern(5012),
	format_hash_streams(5012),
	filter_hash_streams(5012),
	buffer_ppa(100000),
	buffer_general(300000),
	encoder(nullptr),
	vcf_header(nullptr),
	info_containers(new stream_container[100]),
	format_containers(new stream_container[100])
{
	for(U32 i = 0; i < 100; ++i){
		this->info_containers[i].resize(65536*4);
		this->format_containers[i].resize(65536*1000);
	}
}

TomahawkImportWriter::~TomahawkImportWriter(){
	delete this->encoder;
	this->buffer_encode_rle.deleteAll();
	this->buffer_meta.deleteAll();
	this->buffer_encode_simple.deleteAll();
	this->buffer_metaComplex.deleteAll();
	this->buffer_general.deleteAll();
}

bool TomahawkImportWriter::Open(const std::string output){
	this->filename = output;
	this->CheckOutputNames(output);
	this->streamTomahawk.open(this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX, std::ios::out | std::ios::binary);
	this->streamTotempole.open(this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX + '.' + Constants::OUTPUT_INDEX_SUFFIX, std::ios::out | std::ios::binary);

	// Check streams
	if(!this->streamTomahawk.good()){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX << "!" << std::endl;
		return false;
	}
	if(!this->streamTotempole.good()){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX + '.' + Constants::OUTPUT_INDEX_SUFFIX << "!" << std::endl;
		return false;
	}

	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX << "..." << std::endl;
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX + '.' + Constants::OUTPUT_INDEX_SUFFIX << "..." << std::endl;
	}

	// Write Tomahawk and Totempole headers
	this->WriteHeaders();

	// Determine flush limit
	this->DetermineFlushLimit();

	return true;
}

void TomahawkImportWriter::DetermineFlushLimit(void){
	this->flush_limit = this->vcf_header->samples * this->n_variants_limit / 10; // Worst case
	if(this->vcf_header->samples <= Constants::UPPER_LIMIT_SAMPLES_8B - 1)
		this->flush_limit *= sizeof(BYTE);
	else if(this->vcf_header->samples <= Constants::UPPER_LIMIT_SAMPLES_16B - 1)
		this->flush_limit *= sizeof(U16);
	else if(this->vcf_header->samples <= Constants::UPPER_LIMIT_SAMPLES_32B - 1)
		this->flush_limit *= sizeof(U32);
	else this->flush_limit *= sizeof(U64);
}

bool TomahawkImportWriter::OpenExtend(const std::string output){
	this->filename = output;
	this->CheckOutputNames(output);
	this->streamTomahawk.open(output, std::ios::in | std::ios::out | std::ios::binary | std::ios::ate);
	this->streamTotempole.open(output + '.' + Constants::OUTPUT_INDEX_SUFFIX, std::ios::in | std::ios::out | std::ios::binary | std::ios::ate);

	// Check streams
	if(!this->streamTomahawk.good()){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not open: " << output << "!" << std::endl;
		return false;
	}
	if(!this->streamTotempole.good()){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not open: " << output + '.' + Constants::OUTPUT_INDEX_SUFFIX << "!" << std::endl;
		return false;
	}

	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Extending: " << output << "..." << std::endl;
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Extending: " << output + '.' + Constants::OUTPUT_INDEX_SUFFIX << "..." << std::endl;
	}

	U64 tempsize = this->streamTomahawk.tellp();
	this->streamTomahawk.seekp(tempsize - sizeof(U64) * Tomahawk::Constants::eof_length);
	tempsize = this->streamTotempole.tellp();
	this->streamTotempole.seekp(tempsize - sizeof(U64) * Tomahawk::Constants::eof_length);

	// Determine flush limit
	this->DetermineFlushLimit();

	return true;
}

void TomahawkImportWriter::WriteHeaders(void){
	if(this->vcf_header == nullptr){
		std::cerr << Helpers::timestamp("ERROR", "INTERNAL") << "Header not set!" << std::endl;
		exit(1);
	}

	this->streamTotempole.write(Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH);
	this->streamTomahawk.write(Constants::WRITE_HEADER_MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH);

	const U64& samples = this->vcf_header->samples;
	Totempole::TotempoleHeader h(samples);
	this->streamTotempole << h;
	Totempole::TotempoleHeaderBase* hB = reinterpret_cast<Totempole::TotempoleHeaderBase*>(&h);
	this->streamTomahawk << *hB;

	// Write the number of contigs
	const U32 n_contigs = this->vcf_header->contigs.size();
	this->streamTotempole.write(reinterpret_cast<const char*>(&n_contigs), sizeof(U32));

	// Write contig data to Totempole
	// length | n_char | chars[0 .. n_char - 1]
	for(U32 i = 0; i < this->vcf_header->contigs.size(); ++i){
		Totempole::TotempoleContigBase contig(this->vcf_header->contigs[i].length,
											  this->vcf_header->contigs[i].name.size(),
											  this->vcf_header->contigs[i].name);

		this->streamTotempole << contig;
	}

	// Write sample names
	// n_char | chars[0..n_char - 1]
	for(U32 i = 0; i < samples; ++i){
		const U32 n_char = this->vcf_header->sampleNames[i].size();
		this->streamTotempole.write(reinterpret_cast<const char*>(&n_char), sizeof(U32));
		this->streamTotempole.write(reinterpret_cast<const char*>(&this->vcf_header->sampleNames[i][0]), n_char);
	}

	// Push in VCF header and executed line
	buffer_type temp(this->vcf_header->literal_lines.size()*65536);
	for(U32 i = 0; i < this->vcf_header->literal_lines.size(); ++i)
		temp += this->vcf_header->literal_lines[i] + '\n';

	const std::string command = "##tomahawk_importCommand=" + std::string(Constants::LITERAL_COMMAND_LINE)
		+ "; VERSION=" + std::string(VERSION)
		+ "; Date=" + Tomahawk::Helpers::datetime() + "; SIMD=" + SIMD_MAPPING[SIMD_VERSION];

	temp += command;
	this->gzip_controller.Deflate(temp);
	this->streamTotempole.write(&this->gzip_controller.buffer.data[0], this->gzip_controller.buffer.pointer);
	this->gzip_controller.Clear();
	temp.deleteAll();
}

void TomahawkImportWriter::WriteFinal(void){
	// Write EOF
	for(U32 i = 0; i < Constants::eof_length; ++i){
		this->streamTotempole.write(reinterpret_cast<const char*>(&Constants::eof[i]), sizeof(U64));
		this->streamTomahawk.write(reinterpret_cast<const char*>(&Constants::eof[i]), sizeof(U64));
	}

	// Re-open file and overwrite block counts and offset
	const U32 shift = Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH + sizeof(float) + sizeof(U64) + sizeof(BYTE);
	this->streamTotempole.flush();
	std::fstream streamTemp(this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX + '.' + Constants::OUTPUT_INDEX_SUFFIX, std::ios_base::binary | std::ios_base::out | std::ios_base::in);

	if(!streamTemp.good()){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not re-open file!" << std::endl;
		exit(1);
	}

	streamTemp.seekg(shift);
	streamTemp.write(reinterpret_cast<const char*>(&this->n_blocksWritten), sizeof(U32));
	streamTemp.write(reinterpret_cast<const char*>(&this->largest_uncompressed_block), sizeof(U32));
	streamTemp.flush();
	streamTemp.close();
}

void TomahawkImportWriter::setHeader(VCF::VCFHeader& header){
	this->vcf_header = &header;
	this->encoder = new encoder_type(header.samples);
	//this->encoder->DetermineBitWidth();
}

bool TomahawkImportWriter::add(const VCF::VCFLine& line){
	const U32 meta_start_pos = this->buffer_meta.pointer;
	const U32 rle_start_pos  = this->buffer_encode_rle.pointer;
	if(!this->encoder->Encode(line, this->buffer_meta, this->buffer_encode_rle)){
		this->buffer_meta.pointer = meta_start_pos; // reroll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // reroll back
		return false;
	}

	//const U64 n_runs = (this->buffer_encode_rle.pointer - rle_start_pos)/this->rleController_->getBitWidth();
	const meta_base_type& base_meta = *reinterpret_cast<const meta_base_type* const>(&this->buffer_meta[meta_start_pos]);

	if(base_meta.HWE_P < this->filter.HWE_P){
		this->buffer_meta.pointer = meta_start_pos; // reroll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // reroll back
		//std::cerr << "HWE_P < " << this->filter.HWE_P << ": " << base_meta.HWE_P << '\t' << base_meta << std::endl;
		return false;
	}

	if(base_meta.AF < this->filter.MGF){
		this->buffer_meta.pointer = meta_start_pos; // reroll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // reroll back
		//std::cerr << "MAF < " << this->filter.MGF << ": " << base_meta.AF << '\t' << base_meta << std::endl;
		return false;
	}

	if(this->totempole_entry.minPosition == 0)
		this->totempole_entry.minPosition = line.position;

	this->totempole_entry.maxPosition = line.position;
	++this->totempole_entry.n_variants;

	return true;
}

bool TomahawkImportWriter::add(bcf_entry_type& entry, const U32* const ppa){
	// Keep positions
	// If the entry needs to be filtered out
	// then we roll back to these positions
	// In practice we simply move the pointer back
	const U64 meta_start_pos = this->buffer_meta.pointer;
	const U64 simple_start_pos = this->buffer_encode_simple.pointer;
	const U64 rle_start_pos  = this->buffer_encode_rle.pointer;
	meta_base_type meta;

	// Perform run-length encoding
	U64 n_runs = 0;
	if(!this->encoder->Encode(entry, meta, this->buffer_encode_rle, this->buffer_encode_simple, n_runs, ppa)){
		this->buffer_meta.pointer = meta_start_pos; // roll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // roll back
		this->buffer_encode_simple.pointer = simple_start_pos;
		return false;
	}

	// If filtered out: reset pointers
	// Hardy-Weinberg P filter
	if(meta.HWE_P < this->filter.HWE_P){
		this->buffer_meta.pointer = meta_start_pos; // roll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // roll back
		this->buffer_encode_simple.pointer = simple_start_pos;
		//std::cerr << "HWE_P < " << this->filter.HWE_P << ": " << meta.HWE_P << std::endl;
		return false;
	}

	// Minor-genotype frequency filter
	if(meta.AF < this->filter.MGF){
		this->buffer_meta.pointer = meta_start_pos; // roll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // roll back
		this->buffer_encode_simple.pointer = simple_start_pos;
		//std::cerr << "MAF < " << this->filter.MGF << ": " << meta.MGF << std::endl;
		return false;
	}

	// Parse BCF
	this->parseBCF(meta, entry);

	// If the current minPosition is 0
	// then this is the first entry we've seen
	// in this contig. Keep the current position
	// as the last one we've seen
	if(this->totempole_entry.minPosition == 0)
		this->totempole_entry.minPosition = entry.body->POS + 1;

	// Update max position
	this->totempole_entry.maxPosition = entry.body->POS + 1;

	// Push meta to buffer
	// update complex offset position
	meta.virtual_offset_cold_meta = this->buffer_metaComplex.pointer;
	this->buffer_meta += meta;

	// RLE using this word size
	U32 w = ceil(ceil(log2(this->vcf_header->samples + 1))/8);
	if((w > 2) & (w < 4)) w = 4;
	else if(w > 4) w = 8;

	switch(w){
	case 1: this->buffer_meta += (BYTE)n_runs; break;
	case 2: this->buffer_meta += (U16)n_runs; break;
	case 4: this->buffer_meta += (U32)n_runs; break;
	case 8: this->buffer_meta += (U64)n_runs; break;
	default:
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Illegal word-size!" << std::endl;
		exit(1); // unrecoverable error
	}

	// Complex meta data
	Support::TomahawkSupport test;
	if(!test.write(entry, this->buffer_metaComplex)){
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Failed to write complex meta!" << std::endl;
		return false;
	}

	// Update number of entries in block
	++this->totempole_entry.n_variants;

	return true;
}

bool TomahawkImportWriter::parseBCF(meta_base_type& meta, bcf_entry_type& entry){
	//std::cerr << this->body->CHROM << ':' << this->body->POS+1 << '\t' << this->body->n_allele << '\t' << this->body->n_fmt << '\t' << this->body->n_info << '\t' << this->body->n_sample << std::endl;
	U32 internal_pos = entry.filter_start;

	// At FILTER
	// Typed vector
	const bcf_entry_type::base_type& filter_key = *reinterpret_cast<const bcf_entry_type::base_type* const>(&entry.data[internal_pos++]);
	U32 n_filter = filter_key.high;
	if(n_filter == 15) n_filter = entry.getInteger(filter_key.low, internal_pos);
	entry.n_filter = n_filter;
	entry.filter_key = filter_key;

	S32 val = 0;
	while(entry.nextFilter(val, internal_pos)){
		// Hash FILTER value
		// Filter fields have no values
		U32* hash_map_ret = nullptr;
		U32 temp = val;
		if(this->filter_hash_streams.GetItem(&temp, hash_map_ret, sizeof(U32))){
			// exists
		} else {
			this->filter_hash_streams.SetItem(&temp, this->filter_hash_value_counter, sizeof(U32));
			this->filter_values.push_back(val);
			++this->filter_hash_value_counter;
		}
	}

	// At INFO
	U32 info_length;
	BYTE info_value_type;
	while(entry.nextInfo(val, info_length, info_value_type, internal_pos)){
		// Hash INFO values
		U32* hash_map_ret = nullptr;
		U32 mapID = 0;
		U32 temp = val;
		if(this->info_hash_streams.GetItem(&temp, hash_map_ret, sizeof(U32))){
			mapID = *hash_map_ret;
		} else {
			this->info_hash_streams.SetItem(&temp, this->info_hash_value_counter, sizeof(U32));
			mapID = this->info_hash_value_counter;
			this->info_values.push_back(val);
			++this->info_hash_value_counter;
		}

		//
		stream_container& target_container = this->info_containers[mapID];
		if(this->info_containers[mapID].n_entries == 0){
			target_container.setStrideSize(info_length);
			// Set all integer types to U32
			// Change to smaller type later if required
			if(info_value_type == 0)      target_container.setType(4);
			else if(info_value_type == 1) target_container.setType(4);
			else if(info_value_type == 2) target_container.setType(4);
			else if(info_value_type == 3) target_container.setType(4);
			else if(info_value_type == 5) target_container.setType(7);
			else if(info_value_type == 7) target_container.setType(0);
		}
		++target_container;
		if(!target_container.checkStrideSize(info_length))
			target_container.setMixedStrides();

		target_container.addStride(info_length);

		// Flags and integers
		if(info_value_type <= 3){
			for(U32 j = 0; j < info_length; ++j){
				target_container += entry.getInteger(info_value_type, internal_pos);
			}
		}
		// Floats
		else if(info_value_type == 5){
			for(U32 j = 0; j < info_length; ++j){
				target_container += entry.getFloat(internal_pos);
			}
		}
		// Chars
		else if(info_value_type == 7){
			for(U32 j = 0; j < info_length; ++j){
				target_container += entry.getChar(internal_pos);
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible: " << (int)info_value_type << std::endl;
			exit(1);
		}
	}

#if BCF_ASSERT == 1
	// Assert all FILTER and INFO data have been successfully
	// parsed. This is true when the byte pointer equals the
	// start position of the FORMAT fields which are encoded
	// in the meta header structure
	assert(internal_pos == (entry.body->l_shared + sizeof(U32)*2));
#endif

	while(entry.nextFormat(val, info_length, info_value_type, internal_pos)){
		// Hash FORMAT values
		U32* hash_map_ret = nullptr;
		U32 mapID = 0;
		U32 temp = val;
		if(this->format_hash_streams.GetItem(&temp, hash_map_ret, sizeof(U32))){
			mapID = *hash_map_ret;
		} else {
			this->format_hash_streams.SetItem(&temp, this->format_hash_value_counter, sizeof(U32));
			mapID = this->format_hash_value_counter;
			this->format_values.push_back(val);
			std::cerr << Helpers::timestamp("DEBUG") << val << '\t' << info_length << '\t' << (U32)info_value_type << std::endl;
			++this->format_hash_value_counter;
		}

		if(this->format_containers[mapID].n_entries == 0){
			this->format_containers[mapID].setStrideSize(info_length);
		}

		if(mapID == 0){
			switch(info_value_type){
			case 1: internal_pos += this->vcf_header->samples * sizeof(SBYTE) * info_length; break;
			case 2: internal_pos += this->vcf_header->samples * sizeof(S16)   * info_length; break;
			case 3: internal_pos += this->vcf_header->samples * sizeof(S32)   * info_length; break;
			}
			continue;
		}

		stream_container& target_container = this->format_containers[mapID];
		if(this->format_containers[mapID].n_entries == 0){
			target_container.setStrideSize(info_length);
			// Set all integer types to U32
			// Change to smaller type later if required
			if(info_value_type == 0)      target_container.setType(4);
			else if(info_value_type == 1) target_container.setType(4);
			else if(info_value_type == 2) target_container.setType(4);
			else if(info_value_type == 3) target_container.setType(4);
			else if(info_value_type == 5) target_container.setType(7);
			else if(info_value_type == 7) target_container.setType(0);
		}
		++target_container;
		if(!target_container.checkStrideSize(info_length))
			target_container.setMixedStrides();

		target_container.addStride(info_length);

		// Flags and integers
		if(info_value_type <= 3){
			for(U32 j = 0; j < this->vcf_header->samples*info_length; ++j){
				target_container += entry.getInteger(info_value_type, internal_pos);
			}
		}
		// Floats
		else if(info_value_type == 5){
			for(U32 j = 0; j < this->vcf_header->samples*info_length; ++j){
				target_container += entry.getFloat(internal_pos);
			}
		}
		// Chars
		else if(info_value_type == 7){
			for(U32 j = 0; j < this->vcf_header->samples*info_length; ++j){
				target_container += entry.getChar(internal_pos);
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible: " << (int)info_value_type << std::endl;
			exit(1);
		}
	}

	// Hash FILTER pattern
	U32 mapID = 0;
	U32* hash_map_ret = nullptr;
	const U64 hash_filter_vector = entry.hashFilter();
	if(this->filter_hash_pattern.GetItem(&hash_filter_vector, hash_map_ret, sizeof(U64))){
		mapID = *hash_map_ret;
	} else {
		this->filter_hash_pattern.SetItem(&hash_filter_vector, this->filter_hash_pattern_counter, sizeof(U64));
		filter_patterns.push_back(std::vector<U32>());
		for(U32 i = 0; i < entry.filterPointer; ++i){
			this->filter_patterns[this->filter_hash_pattern_counter].push_back(entry.filterID[i]);
		}
		assert(this->filter_hash_pattern_counter < 65536);
		mapID = this->filter_hash_pattern_counter;
		++this->filter_hash_pattern_counter;
	}
	// Store this map in the meta
	meta.FILTER_map_ID = mapID;

	// Hash INFO pattern
	mapID = 0;
	hash_map_ret = nullptr;
	// Hash INFO vector of identifiers
	const U64 hash_info_vector = entry.hashInfo();
	if(this->info_hash_pattern.GetItem(&hash_info_vector, hash_map_ret, sizeof(U64))){
		mapID = *hash_map_ret;
	} else {
		this->info_hash_pattern.SetItem(&hash_info_vector, this->info_hash_pattern_counter, sizeof(U64));
		this->info_patterns.push_back(std::vector<U32>());
		for(U32 i = 0; i < entry.infoPointer; ++i){
			this->info_patterns[this->info_hash_pattern_counter].push_back(entry.infoID[i]);
		}
		assert(this->info_hash_pattern_counter < 65536);
		mapID = this->info_hash_pattern_counter;
		++this->info_hash_pattern_counter;
	}
	// Store this map in the meta
	meta.INFO_map_ID = mapID;

	// Hash FORMAT pattern
	mapID = 0;
	hash_map_ret = nullptr;
	const U64 hash_format_vector = entry.hashFormat();
	if(this->format_hash_pattern.GetItem(&hash_format_vector, hash_map_ret, sizeof(U64))){
		mapID = *hash_map_ret;
	} else {
		this->format_hash_pattern.SetItem(&hash_format_vector, this->format_hash_pattern_counter, sizeof(U64));
		format_patterns.push_back(std::vector<U32>());
		for(U32 i = 0; i < entry.formatPointer; ++i){
			this->format_patterns[this->format_hash_pattern_counter].push_back(entry.formatID[i]);
		}
		assert(this->format_hash_pattern_counter < 65536);
		mapID = this->format_hash_pattern_counter;
		++this->format_hash_pattern_counter;
	}
	// Store this map in the meta
	meta.FORMAT_map_ID = mapID;

	// Return
	return true;
}

// flush and write
bool TomahawkImportWriter::flush(Algorithm::RadixSortGT& permuter){
	if(this->buffer_meta.size() == 0)
		return false;

	// Update sizes of streams
	this->totempole_entry.byte_offset = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
	this->totempole_entry.l_meta = this->buffer_meta.pointer;
	this->totempole_entry.l_meta_complex = this->buffer_metaComplex.pointer;
	this->totempole_entry.l_gt_rle = this->buffer_encode_rle.pointer;
	this->totempole_entry.l_gt_simple = this->buffer_encode_simple.pointer;

	// Todo:
	// Build map
	// Construct bitvectors
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "INFO: " << this->info_values.size() << " entries -> " << ceil((float)this->info_values.size()/8) << " bytes" << std::endl;
	if(!this->constructBitVector(this->info_hash_streams, this->info_values, this->info_patterns)){
		std::cerr << "failed bit" << std::endl;
		return false;
	}

	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "FORMAT: " << this->format_values.size() << " entries -> " << ceil((float)this->format_values.size()/8) << " bytes" << std::endl;
	if(!this->constructBitVector(this->format_hash_streams, this->format_values, this->format_patterns)){
		std::cerr << "failed bit" << std::endl;
		return false;
	}

	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "FILTER: " << this->filter_values.size() << " entries -> " << ceil((float)this->filter_values.size()/8) << " bytes" << std::endl;
	if(!this->constructBitVector(this->filter_hash_streams, this->filter_values, this->filter_patterns)){
		std::cerr << "failed bit" << std::endl;
		return false;
	}

	//std::cerr << this->totempole_entry.n_variants << '\t' << this->buffer_meta.pointer << '\t' << this->buffer_metaComplex.pointer << '\t' << this->buffer_encode_rle.pointer << '\t' << this->buffer_encode_simple.pointer << std::endl;

	// Merge GT data
	/*
	U32 freq[256];
	memset(&freq[0], 0, sizeof(U32)*256);
	for(U32 i = 0; i < this->buffer_encode_rle.size(); ++i){
		++freq[(BYTE)this->buffer_encode_rle[i]];
	}

	double p_sum = 0;
	for(U32 i = 0; i < 256; ++i){
		double p = (double)freq[i]/this->buffer_encode_rle.size();
		if(p == 0) continue;

		std::cout << i << '\t' << freq[i] << '\t' << -p*log2(p) << std::endl;
		p_sum += p*log2(p);
	}
	*/

	// Test to iterate over complex

	//U32 p = 0;
	/*
	TomahawkBlockIterator<U16> b(this->buffer_meta.data, this->buffer_meta.pointer, this->totempole_entry);
	while(true){
		//std::cerr << b.getMeta() << std::endl;
		//std::cerr << p++ << '/' << this->totempole_entry.n_variants << ": " << b.isRLE() << '\t' << b.size() << '\t' << b.getMeta() << std::endl;

		if(b.getMeta().controller.biallelic == 0 && b.getMeta().controller.rle == 0){
			std::cerr << b.getMeta().position << '\t' << b.getMeta().n_runs << std::endl;
			const BYTE n_alleles = b.getMetaComplex().n_allele + 1;
			if(n_alleles < 8){
				const Support::TomahawkRunSimple<BYTE>* field = nullptr;
				while(b.nextRunSimple(field)){
					std::cerr << *field << '\t';
				}
				std::cerr << std::endl;

				//exit(1);
			} else {
				std::cerr << this->totempole_entry.n_variants << ": " << b.isRLE() << '\t' << b.size() << '\t' << b.getMeta() << std::endl;
				const Support::TomahawkRunSimple<U16>* field = nullptr;
				while(b.nextRunSimple(field)){
					//std::cerr << *reinterpret_cast<const U16*>(field) << '\t' << *field << std::endl;
					std::cerr << *field << '\t';
				}
				std::cerr << std::endl;
				//exit(1);
			}
		} else if(b.getMeta().controller.biallelic == 0 && b.getMeta().controller.rle == 1){
			const BYTE& rle_type = b.getMeta().controller.rle_type;
			const BYTE shift_size = ceil(log2(b.getMetaComplex().n_allele + 1));
			const BYTE allele_mask = (1 << shift_size) - 1;

			if(rle_type == 0){
				std::cerr << b.getMeta().position << std::endl;
				std::cerr << "nallelic and rle: " << (int)rle_type << " n_alleles: " << b.getMetaComplex().n_allele << std::endl;
				std::cerr << (int)shift_size << ',' << (int)rle_type << " mask: " << std::bitset<8>(allele_mask) << std::endl;

				const BYTE* field = nullptr;

				while(b.nextRunSimple(field)){
					const U32 l = ((*field >> (1 + 2*shift_size)));
					for(U32 i = 0; i < l; ++i)
						std::cerr << (int)(((*field >> 1) & allele_mask) - 1) << ((*field&1) ? '|' : '/') << (int)(((*field >> (1 + shift_size)) & allele_mask) - 1) << '\t';
					//std::cerr << "A: " << (int)((*field >> 1) & allele_mask) << '\t'
					//		  << "B: " << (int)((*field >> (1 + shift_size)) & allele_mask) << '\t'
					//		  << "run; " << (int)((*field >> (1 + 2*shift_size))) << '\t'
					//		  << std::bitset<8>(*field) << std::endl;
				}
				exit(1);
			}
		}


		else if(b.getMeta().controller.biallelic == 1 && b.getMeta().controller.rle == 1) {
			/*
			const BYTE& rle_type = b.getMeta().controller.rle_type;
			const BYTE shift_size = 1 + b.getMeta().controller.anyMissing;
			const BYTE mixed_phasing = b.getMeta().controller.mixed_phasing;
			const BYTE allele_mask = (1 << shift_size) - 1;

			if(rle_type == 1){
				const U16* field = nullptr;
				std::cerr << (int)shift_size << ',' << (int)mixed_phasing << ',' << (int)rle_type << " mask: " << std::bitset<16>(allele_mask) << std::endl;
				while(b.nextRun(field)){
					std::cerr << "A: " << (int)((*field >> mixed_phasing) & allele_mask) << '\t'
							  << "B: " << (int)((*field >> (mixed_phasing + shift_size)) & allele_mask) << '\t'
							  << "run; " << (int)((*field >> (mixed_phasing + 2*shift_size))) << '\t'
							  << std::bitset<16>(*field) << std::endl;
				}
				exit(1);
			}
			*/
		/*}


		if(!++b) break;
	}
	*/

	// Split U32 values into 4 streams
	const U32 partition = this->vcf_header->samples;
	//buffer_type test(partition*sizeof(U32)*10);
	bytePreprocessor(permuter.getPPA(), partition, this->buffer_general.data);

	// test bitshuffle
	memset(this->buffer_ppa.data, 0, 2*sizeof(BYTE)*partition);
	bytePreprocessBits(&this->buffer_general.data[partition*2], partition, this->buffer_ppa.data);
	bytePreprocessBits(&this->buffer_general.data[partition*3], partition, &this->buffer_ppa.data[partition]);
	//const BYTE* const test_bit = reinterpret_cast<const BYTE* const>(this->buffer_ppa.data);
	//for(U32 i = 0; i < partition; ++i){
//		std::cerr << (U32)test_bit[i] << ' ';
	//}
	//std::cerr << std::endl;
	this->buffer_ppa.pointer = 2*partition;
	this->gzip_controller.Deflate(this->buffer_ppa);
	this->streamTomahawk << this->gzip_controller;
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "PPA\t" << 2*partition << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)2*partition/this->gzip_controller.buffer.size() << std::endl;
	std::cout << partition*2 << '\t' << this->gzip_controller.buffer.size() << '\t';
	this->gzip_controller.Clear();
	this->buffer_general.reset();

	/*

	this->buffer_ppa.resize(partition*8);
	U32 ret_sep = 0;
	buffer_type temp2(partition*8);
	buffer_type temp3(partition*8);
	std::cerr << "size: " << temp2.capacity() << std::endl;
	for(U32 i = 0; i < 4; ++i){
		std::cerr << "compressing: " << partition*i << "->" << partition*(i+1) << " at " << this->buffer_ppa.pointer << '/' << this->buffer_ppa.capacity() << std::endl;
		size_t ret = FSE_compress(&this->buffer_ppa.data[this->buffer_ppa.pointer], this->buffer_ppa.capacity(),
				                  &test.data[partition*i], partition);
		std::cerr << "compressed: " << ret << " input: " << partition << std::endl;
		if(FSE_isError(ret) || ret == 0){
			std::cerr << "failed to compress: " << ret << std::endl;
			exit(1);
		}

		if(ret == 1)
			continue;

		size_t retD = FSE_decompress(temp2.data, temp2.capacity(), &this->buffer_ppa.data[this->buffer_ppa.pointer], ret);
		if(FSE_isError(retD)){
			std::cerr << "decompress error: "  << retD << std::endl;
			exit(1);
		}
		std::cerr << "decompressed: " << retD << std::endl;
		ret_sep += ret;
		this->buffer_ppa.pointer += ret;
	}
	this->streamTomahawk << this->buffer_ppa;
	std::cerr << Helpers::timestamp("DEBUG","IMPORT") << "PPA\t" << this->vcf_header->samples*sizeof(U32) << '\t' << this->buffer_ppa.pointer << std::endl;
	this->buffer_ppa.reset();
	*/

	// Merge data to single buffer
	// and compress
	//
	// Meta
	//this->buffer_meta += this->buffer_metaComplex;
	this->gzip_controller.Deflate(this->buffer_meta); // Deflate block
	this->streamTomahawk << this->gzip_controller; // Write tomahawk output
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "META\t" << this->buffer_meta.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->buffer_meta.size() / this->gzip_controller.buffer.size() << std::endl;
	std::cout << this->buffer_meta.size() << '\t' << this->gzip_controller.buffer.size() << '\t';
	this->gzip_controller.Clear(); // Clean up gzip controller

	this->gzip_controller.Deflate(this->buffer_metaComplex); // Deflate block
	this->streamTomahawk << this->gzip_controller; // Write tomahawk output
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "META-C\t" << this->buffer_metaComplex.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->buffer_metaComplex.size() / this->gzip_controller.buffer.size() << std::endl;
	std::cout << this->buffer_metaComplex.size() << '\t' << this->gzip_controller.buffer.size() << '\t';
	this->gzip_controller.Clear(); // Clean up gzip controller

	// RLE
	//this->buffer_encode_rle += this->buffer_encode_simple;
	this->gzip_controller.Deflate(this->buffer_encode_rle); // Deflate block
	this->streamTomahawk << this->gzip_controller; // Write tomahawk output
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "RLE\t" << this->buffer_encode_rle.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->buffer_encode_rle.size() / this->gzip_controller.buffer.size() << std::endl;
	std::cout << this->buffer_encode_rle.size() << '\t' << this->gzip_controller.buffer.size() << '\t';
	this->gzip_controller.Clear(); // Clean up gzip controller

	// If there is any simple data
	if(this->buffer_encode_simple.size() > 0){
		this->gzip_controller.Deflate(this->buffer_encode_simple); // Deflate block
		this->streamTomahawk << this->gzip_controller; // Write tomahawk output
		std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "RLE-S\t" << this->buffer_encode_simple.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->buffer_encode_simple.size() / this->gzip_controller.buffer.size() << std::endl;
		std::cout << this->buffer_encode_simple.size() << '\t' << this->gzip_controller.buffer.size() << '\t';
		this->gzip_controller.Clear(); // Clean up gzip controller
	}

	// Keep track of largest block observed
	if(this->buffer_encode_rle.size() > this->largest_uncompressed_block)
		this->largest_uncompressed_block = this->buffer_encode_rle.size();

	//this->totempole_entry.l_uncompressed = this->buffer_meta.size(); // Store uncompressed size
	this->totempole_entry.byte_offset_end = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
	this->streamTotempole << this->totempole_entry; // Write totempole output
	++this->n_blocksWritten; // update number of blocks written
	this->n_variants_written += this->totempole_entry.n_variants; // update number of variants written
	this->totempole_entry.reset();

	// Dispatch values into streams
	for(U32 i = 0; i < this->info_hash_value_counter; ++i){
		S32 ret_size = this->recodeStream(this->info_containers[i]);
		if(ret_size == -1){
			std::cerr << "failed recode @ " << i << std::endl;
			exit(1);
		}
		std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "INFO " << this->info_values[i] << '\t' << ret_size << '\t' << (float)this->info_containers[i].buffer_data.size()/ret_size << '\t' << this->info_containers[i].addSize << std::endl;
		if(this->info_containers[i].addSize == -1){
			this->gzip_controller.Deflate(this->info_containers[i].buffer_strides);
			std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "INFO-ADD " << this->info_containers[i].buffer_strides.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->info_containers[i].buffer_strides.size()/this->gzip_controller.buffer.size() << std::endl;
			this->streamTomahawk << this->gzip_controller;
			this->gzip_controller.Clear();
		}
	}

	// Dispatch values into streams
	for(U32 i = 0; i < this->format_hash_value_counter; ++i){
		S32 ret_size = this->recodeStream(this->format_containers[i]);
		if(ret_size == -1){
			std::cerr << "failed recode @ " << i << std::endl;
			exit(1);
		}
		std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "FORMAT " << this->format_values[i] << '\t' << ret_size << '\t' << (float)this->format_containers[i].buffer_data.size()/ret_size << '\t' << this->format_containers[i].addSize << std::endl;
		if(this->format_containers[i].addSize == -1){
			this->gzip_controller.Deflate(this->format_containers[i].buffer_strides);
			std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "FORMAT-ADD " << this->format_containers[i].buffer_strides.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->format_containers[i].buffer_strides.size()/this->gzip_controller.buffer.size() << std::endl;
			this->streamTomahawk << this->gzip_controller;
			this->gzip_controller.Clear();
		}
	}

	this->reset(); // reset buffers

	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "END FLUSH" << std::endl;

	return true;
}

bool TomahawkImportWriter::constructBitVector(hash_table& htable, const id_vector& values, const pattern_vector& patterns){
	const BYTE bitvector_width = ceil((float)values.size()/8);
	BYTE* bitvector = new BYTE[bitvector_width];
	memset(bitvector, 0, sizeof(BYTE)*bitvector_width);

	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "Mapping..." << std::endl;

	// Cycle over pattern size
	for(U32 i = 0; i < patterns.size(); ++i){
		// Dump data
		std::cerr << i << '\t';
		for(U32 j = 0; j < patterns[i].size(); ++j){
			std::cerr << patterns[i][j] << '\t';
		}
		std::cerr << std::endl;

		//
		std::cerr << i << '\t';
		for(U32 j = 0; j < patterns[i].size(); ++j){
			U32* retval = nullptr;
			if(!htable.GetItem(&patterns[i][j], retval, sizeof(U32))){
				std::cerr << "impossible" << std::endl;
				exit(1);
			}
			bitvector[*retval/8] ^= 1 << (*retval % 8);
			std::cerr << *retval << '\t';
		}
		std::cerr << std::endl;

		std::cerr << i << '\t';
		for(U32 j = 0; j < bitvector_width; ++j)
			std::cerr << std::bitset<8>(bitvector[j]);

		std::cerr << std::endl << std::endl;
		memset(bitvector, 0, sizeof(BYTE)*bitvector_width);
	}
	std::cerr << std::endl;
	delete [] bitvector;

	return true;
}

S32 TomahawkImportWriter::recodeStream(stream_container& stream){
	S32 ret_size = -1;

	if(stream.stream_data_type == 4){
		const S32* dat = reinterpret_cast<const S32*>(stream.buffer_data.data);
		S32 min = dat[0];
		S32 max = dat[0];
		S32 prev_value = dat[0];
		bool is_uniform = true;

		for(U32 j = 1; j < stream.n_entries; ++j){
			if(dat[j] < min) min = dat[j];
			if(dat[j] > max) max = dat[j];
			if(prev_value != dat[j]) is_uniform = false;
			prev_value = dat[j];
		}

		BYTE byte_width = 0;
		if(min < 0) byte_width = ceil((ceil(log2(abs(min) + 1))+1)/8);  // One bit is used for sign
		else byte_width = ceil(ceil(log2(max + 1))/8);

		if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
		else if(byte_width > 4) byte_width = 8;
		if(byte_width == 0) byte_width = 1;

		// Phase 2
		// Here we re-encode values using the smallest possible
		// word-size
		if(is_uniform){
			// Non-negative
			if(min >= 0){
				switch(byte_width){
				case 1: this->buffer_general += (BYTE)min; break;
				case 2: this->buffer_general += (U16)min; break;
				case 4: this->buffer_general += (U32)min; break;
				case 8: this->buffer_general += (U64)min; break;
				default: std::cerr << "illegal: " << std::endl; exit(1);
				}
			} else {
				switch(byte_width){
				case 1: this->buffer_general += (SBYTE)min; break;
				case 2: this->buffer_general += (S16)min; break;
				case 4: this->buffer_general += (S32)min; break;
				default: std::cerr << "illegal" << std::endl; exit(1);
				}
			}
		} else {
			//std::cerr << "non-uniform" << std::endl;
			dat = reinterpret_cast<const S32*>(stream.buffer_data.data);
			// Is non-negative
			if(min >= 0){
				if(byte_width == 1){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (BYTE)*(dat++);
				} else if(byte_width == 2){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (U16)*(dat++);
				} else if(byte_width == 4){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (U32)*(dat++);
				} else if(byte_width == 8){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (U64)*(dat++);
				} else {
					std::cerr << "illegal" << std::endl;
					exit(1);
				}
			}
			// Is negative
			else {
				if(byte_width == 1){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (SBYTE)*(dat++);
				} else if(byte_width == 2){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (S16)*(dat++);
				} else if(byte_width == 4){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (S32)*(dat++);
				} else {
					std::cerr << "illegal" << std::endl;
					exit(1);
				}
			}
		}

		min = 0; max = 0;
		is_uniform = true;

		if(is_uniform){
			this->streamTomahawk.write(this->buffer_general.data, byte_width);
			ret_size = byte_width;
		} else {
			this->gzip_controller.Deflate(this->buffer_general);
			this->streamTomahawk << this->gzip_controller;
			ret_size = this->gzip_controller.buffer.size();
			this->gzip_controller.Clear();
		}
		this->buffer_general.reset();

	}
	// Is not an integer
	else {
		this->gzip_controller.Deflate(stream.buffer_data);
		this->streamTomahawk << this->gzip_controller;
		ret_size = this->gzip_controller.buffer.size();
		this->gzip_controller.Clear();
	}

	return(ret_size);
}

void TomahawkImportWriter::CheckOutputNames(const std::string& input){
	std::vector<std::string> paths = Helpers::filePathBaseExtension(input);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == Constants::OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &Constants::OUTPUT_SUFFIX[0], Constants::OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}


} /* namespace Tomahawk */
