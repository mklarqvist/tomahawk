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

namespace Tomahawk {

TomahawkImportWriter::TomahawkImportWriter(const filter_type& filter) :
	flush_limit(1000000),
	n_variants_limit(1024),
	blocksWritten_(0),
	variants_written(0),
	variants_complex_written(0),
	largest_uncompressed_block_(0),
	filter(filter),
	encoder(nullptr),
	buffer_encode_rle(flush_limit*2),
	buffer_encode_simple(flush_limit*2),
	buffer_meta(flush_limit*2),
	buffer_metaComplex(flush_limit*2),
	vcf_header_(nullptr)
{}

TomahawkImportWriter::~TomahawkImportWriter(){
	delete this->encoder;
	this->buffer_encode_rle.deleteAll();
	this->buffer_meta.deleteAll();
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
	this->flush_limit = this->vcf_header_->samples * this->n_variants_limit / 10; // Worst case
	if(this->vcf_header_->samples <= Constants::UPPER_LIMIT_SAMPLES_8B - 1)
		this->flush_limit *= sizeof(BYTE);
	else if(this->vcf_header_->samples <= Constants::UPPER_LIMIT_SAMPLES_16B - 1)
		this->flush_limit *= sizeof(U16);
	else if(this->vcf_header_->samples <= Constants::UPPER_LIMIT_SAMPLES_32B - 1)
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
	this->streamTomahawk.seekp(tempsize - sizeof(U64)*Tomahawk::Constants::eof_length);
	tempsize = this->streamTotempole.tellp();
	this->streamTotempole.seekp(tempsize - sizeof(U64)*Tomahawk::Constants::eof_length);

	// Determine flush limit
	this->DetermineFlushLimit();

	return true;
}

void TomahawkImportWriter::WriteHeaders(void){
	if(this->vcf_header_ == nullptr){
		std::cerr << Helpers::timestamp("ERROR", "INTERNAL") << "Header not set!" << std::endl;
		exit(1);
	}

	this->streamTotempole.write(Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH);
	this->streamTomahawk.write(Constants::WRITE_HEADER_MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH);

	const U64& samples = this->vcf_header_->samples;
	Totempole::TotempoleHeader h(samples);
	this->streamTotempole << h;
	Totempole::TotempoleHeaderBase* hB = reinterpret_cast<Totempole::TotempoleHeaderBase*>(&h);
	this->streamTomahawk << *hB;

	// Write the number of contigs
	const U32 n_contigs = this->vcf_header_->contigs.size();
	this->streamTotempole.write(reinterpret_cast<const char*>(&n_contigs), sizeof(U32));

	// Write contig data to Totempole
	// length | n_char | chars[0 .. n_char - 1]
	for(U32 i = 0; i < this->vcf_header_->contigs.size(); ++i){
		Totempole::TotempoleContigBase contig(this->vcf_header_->contigs[i].length,
											  this->vcf_header_->contigs[i].name.size(),
											  this->vcf_header_->contigs[i].name);

		this->streamTotempole << contig;
	}

	// Write sample names
	// n_char | chars[0..n_char - 1]
	for(U32 i = 0; i < samples; ++i){
		const U32 n_char = this->vcf_header_->sampleNames[i].size();
		this->streamTotempole.write(reinterpret_cast<const char*>(&n_char), sizeof(U32));
		this->streamTotempole.write(reinterpret_cast<const char*>(&this->vcf_header_->sampleNames[i][0]), n_char);
	}

	// Push in VCF header and executed line
	buffer_type temp(this->vcf_header_->literal_lines.size()*65536);
	for(U32 i = 0; i < this->vcf_header_->literal_lines.size(); ++i)
		temp += this->vcf_header_->literal_lines[i] + '\n';

	const std::string command = "##tomahawk_importCommand=" + std::string(Constants::LITERAL_COMMAND_LINE)
		+ "; VERSION=" + std::string(VERSION)
		+ "; Date=" + Tomahawk::Helpers::datetime() + "; SIMD=" + SIMD_MAPPING[SIMD_VERSION];

	temp += command;
	this->gzip_controller_.Deflate(temp);
	this->streamTotempole.write(&this->gzip_controller_.buffer.data[0], this->gzip_controller_.buffer.pointer);
	this->gzip_controller_.Clear();
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
	streamTemp.write(reinterpret_cast<const char*>(&this->blocksWritten_), sizeof(U32));
	streamTemp.write(reinterpret_cast<const char*>(&this->largest_uncompressed_block_), sizeof(U32));
	streamTemp.flush();
	streamTemp.close();
}

void TomahawkImportWriter::setHeader(VCF::VCFHeader& header){
	this->vcf_header_ = &header;
	this->encoder = new encoder_type(header.samples);
	this->encoder->DetermineBitWidth();
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

	if(base_meta.MGF < this->filter.MGF){
		this->buffer_meta.pointer = meta_start_pos; // reroll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // reroll back
		//std::cerr << "MAF < " << this->filter.MAF << ": " << base_meta.MAF << '\t' << base_meta << std::endl;
		return false;
	}

	if(this->totempole_entry.minPosition == 0)
		this->totempole_entry.minPosition = line.position;

	this->totempole_entry.maxPosition = line.position;
	//if(line.c)
	++this->totempole_entry.n_variantsRLE;

	return true;
}

bool TomahawkImportWriter::add(const bcf_entry_type& line){
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
	if(!this->encoder->Encode(line, meta, this->buffer_encode_rle, this->buffer_encode_simple, n_runs)){
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
	if(meta.MGF < this->filter.MGF){
		this->buffer_meta.pointer = meta_start_pos; // roll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // roll back
		this->buffer_encode_simple.pointer = simple_start_pos;
		//std::cerr << "MAF < " << this->filter.MGF << ": " << meta.MGF << std::endl;
		return false;
	}

	// If the current minPosition is 0
	// then this is the first entry we've seen
	// in this contig. Keep the current position
	// as the last one we've seen
	if(this->totempole_entry.minPosition == 0)
		this->totempole_entry.minPosition = line.body->POS + 1;

	// Update max position
	this->totempole_entry.maxPosition = line.body->POS + 1;

	// Push meta to buffer
	// update complex offset position
	meta.virtual_offset_complex = this->buffer_metaComplex.pointer;
	meta.controller.hasComplex = true;
	this->buffer_meta += meta;
	switch(this->encoder->getBitWidth()){
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
	if(!test.write(line, this->buffer_metaComplex)){
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Failed to write complex meta!" << std::endl;
		return false;
	}

	// Update number of entries in block
	if(meta.controller.biallelic)
		++this->totempole_entry.n_variantsRLE;
	else ++this->totempole_entry.n_variantsComplex;

	return true;
}

// flush and write
bool TomahawkImportWriter::flush(void){
	if(this->buffer_meta.size() == 0)
		return false;

	// Update sizes of streams
	this->totempole_entry.l_meta = this->buffer_meta.pointer;
	this->totempole_entry.l_meta_complex = this->buffer_metaComplex.pointer;
	this->totempole_entry.l_rle = this->buffer_encode_rle.pointer;
	this->totempole_entry.l_simple = this->buffer_encode_simple.pointer;

	std::cerr << this->totempole_entry.n_variantsRLE+this->totempole_entry.n_variantsComplex << '\t' << this->totempole_entry.n_variantsRLE << '\t' << this->totempole_entry.n_variantsComplex << '\t' << this->buffer_meta.pointer << '\t' << this->buffer_metaComplex.pointer << '\t' << this->buffer_encode_rle.pointer << '\t' << this->buffer_encode_simple.pointer << '\t';

	this->totempole_entry.byte_offset = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
	//this->buffer_meta += this->buffer_metaComplex;
	//this->buffer_encode_rle += this->buffer_encode_simple;

	this->buffer_meta += this->buffer_encode_rle;
	this->buffer_meta += this->buffer_encode_simple;
	this->buffer_meta += this->buffer_metaComplex;

	this->gzip_controller_.Deflate(this->buffer_meta); // Deflate block
	this->streamTomahawk << this->gzip_controller_; // Write tomahawk output
	std::cerr << this->gzip_controller_.buffer.pointer << std::endl;
	this->gzip_controller_.Clear(); // Clean up gzip controller

	// Keep track of largest block observed
	if(this->buffer_meta.size() > this->largest_uncompressed_block_)
		this->largest_uncompressed_block_ = this->buffer_meta.size();

	this->totempole_entry.l_uncompressed = this->buffer_meta.size(); // Store uncompressed size
	this->totempole_entry.byte_offset_end = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
	this->streamTotempole << this->totempole_entry; // Write totempole output
	++this->blocksWritten_; // update number of blocks written
	this->variants_written += this->totempole_entry.n_variantsRLE + this->totempole_entry.n_variantsComplex; // update number of variants written

	this->reset(); // reset buffers

	return true;
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
