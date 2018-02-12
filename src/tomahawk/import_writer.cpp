#include "import_writer.h"

namespace Tomahawk {

ImportWriter::ImportWriter(const filter_type& filter) :
	flush_limit(1000000),
	n_variants_limit(1024),
	blocksWritten_(0),
	variants_written_(0),
	largest_uncompressed_block_(0),
	filter(filter),
	rleController_(nullptr),
	buffer_rle_(flush_limit*2),
	buffer_meta_(flush_limit*2),
	vcf_header_(nullptr)
{}

ImportWriter::~ImportWriter(){
	delete this->rleController_;
	this->buffer_rle_.deleteAll();
	this->buffer_meta_.deleteAll();
}

bool ImportWriter::Open(const std::string output){
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

void ImportWriter::DetermineFlushLimit(void){
	this->flush_limit = this->vcf_header_->samples * this->n_variants_limit / 10; // Worst case
	if(this->vcf_header_->samples <= Constants::UPPER_LIMIT_SAMPLES_8B - 1)
		this->flush_limit *= sizeof(BYTE);
	else if(this->vcf_header_->samples <= Constants::UPPER_LIMIT_SAMPLES_16B - 1)
		this->flush_limit *= sizeof(U16);
	else if(this->vcf_header_->samples <= Constants::UPPER_LIMIT_SAMPLES_32B - 1)
		this->flush_limit *= sizeof(U32);
	else this->flush_limit *= sizeof(U64);
}

bool ImportWriter::OpenExtend(const std::string output){
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

void ImportWriter::WriteHeaders(void){
	if(this->vcf_header_ == nullptr){
		std::cerr << Helpers::timestamp("ERROR", "INTERNAL") << "Header not set!" << std::endl;
		exit(1);
	}

	this->streamTotempole.write(Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH);
	this->streamTomahawk.write(Constants::WRITE_HEADER_MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH);

	/*
	const U64& samples = this->vcf_header_->samples;
	Totempole::TotempoleHeader h(samples);
	this->streamTotempole << h;
	Totempole::TotempoleHeaderBase* hB = reinterpret_cast<Totempole::TotempoleHeaderBase*>(&h);
	this->streamTomahawk << *hB;
	*/
	std::cerr << "not implemented; fix me" << std::endl;
	exit(1);


	// Write out dummy variable for IO offset
	U32 nothing = 0; // Dummy variable
	size_t posOffset = this->streamTotempole.tellp(); // remember current IO position
	this->streamTotempole.write(reinterpret_cast<const char*>(&nothing), sizeof(U32)); // data offset

	// Write the number of contigs
	const U32 n_contigs = this->vcf_header_->contigs.size();
	this->streamTotempole.write(reinterpret_cast<const char*>(&n_contigs), sizeof(U32));

	// Write contig data to Totempole
	// length | n_char | chars[0 .. n_char - 1]
	for(U32 i = 0; i < this->vcf_header_->contigs.size(); ++i){
		Totempole::HeaderContig contig(this->vcf_header_->contigs[i].length,
											  this->vcf_header_->contigs[i].name.size(),
											  this->vcf_header_->contigs[i].name);

		this->streamTotempole << contig;
	}

	// Write sample names
	// n_char | chars[0..n_char - 1]
	/*
	for(U32 i = 0; i < samples; ++i){
		const U32 n_char = this->vcf_header_->sampleNames[i].size();
		this->streamTotempole.write(reinterpret_cast<const char*>(&n_char), sizeof(U32));
		this->streamTotempole.write(reinterpret_cast<const char*>(&this->vcf_header_->sampleNames[i][0]), n_char);
	}
	*/

	// Push in VCF header and executed line
	buffer_type temp(this->vcf_header_->literal_lines.size()*65536);
	for(U32 i = 0; i < this->vcf_header_->literal_lines.size(); ++i)
		temp += this->vcf_header_->literal_lines[i] + '\n';

	const std::string command = "##tomahawk_importCommand=" + std::string(Constants::LITERAL_COMMAND_LINE)
		+ "; VERSION=" + std::string(VERSION)
		+ "; Date=" + Tomahawk::Helpers::datetime() + "; SIMD=" + SIMD_MAPPING[SIMD_VERSION];

	temp += command;
	this->gzip_controller_.Deflate(temp);
	this->streamTotempole.write(this->gzip_controller_.buffer.data(), this->gzip_controller_.buffer.n_chars);
	this->gzip_controller_.Clear();
	temp.deleteAll();

	U32 curPos = this->streamTotempole.tellp(); // remember current IO position
	this->streamTotempole.seekp(posOffset); // seek to previous position
	this->streamTotempole.write(reinterpret_cast<const char*>(&curPos), sizeof(U32)); // overwrite data offset
	this->streamTotempole.seekp(curPos); // seek back to current IO position
}

void ImportWriter::WriteFinal(void){
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

void ImportWriter::setHeader(VCF::VCFHeader& header){
	this->vcf_header_ = &header;
	this->rleController_ = new Algorithm::GenotypeEncoder(header.samples);
	this->rleController_->DetermineBitWidth();
}

bool ImportWriter::add(const VCF::VCFLine& line){
	const U32 meta_start_pos = this->buffer_meta_.n_chars;
	const U32 rle_start_pos  = this->buffer_rle_.n_chars;
	if(!this->rleController_->RunLengthEncode(line, this->buffer_meta_, this->buffer_rle_)){
		this->buffer_meta_.n_chars = meta_start_pos; // reroll back
		this->buffer_rle_.n_chars  = rle_start_pos; // reroll back
		return false;
	}

	const U64 n_runs = (this->buffer_rle_.n_chars - rle_start_pos)/this->rleController_->getBitWidth();
	const MetaEntryBase& base_meta = *reinterpret_cast<const MetaEntryBase* const>(&this->buffer_meta_[meta_start_pos]);

	if(n_runs == 1){
		this->buffer_meta_.n_chars = meta_start_pos; // reroll back
		this->buffer_rle_.n_chars  = rle_start_pos; // reroll back
		//std::cerr << "singleton" << std::endl;
		return false;
	}

	if(base_meta.HWE_P < this->filter.HWE_P){
		this->buffer_meta_.n_chars = meta_start_pos; // reroll back
		this->buffer_rle_.n_chars  = rle_start_pos; // reroll back
		//std::cerr << "HWE_P < " << this->filter.HWE_P << ": " << base_meta.HWE_P << '\t' << base_meta << std::endl;
		return false;
	}

	if(base_meta.MAF < this->filter.MAF){
		this->buffer_meta_.n_chars = meta_start_pos; // reroll back
		this->buffer_rle_.n_chars  = rle_start_pos; // reroll back
		//std::cerr << "MAF < " << this->filter.MAF << ": " << base_meta.MAF << '\t' << base_meta << std::endl;
		return false;
	}

	if(this->totempole_entry.min_position == 0)
		this->totempole_entry.min_position = line.position;

	this->totempole_entry.max_position = line.position;
	++this->totempole_entry;

	return true;
}

bool ImportWriter::add(const BCF::BCFEntry& line){
	const U32 meta_start_pos = this->buffer_meta_.n_chars;
	const U32 rle_start_pos  = this->buffer_rle_.n_chars;
	if(!this->rleController_->RunLengthEncode(line, this->buffer_meta_, this->buffer_rle_)){
		this->buffer_meta_.n_chars = meta_start_pos; // reroll back
		this->buffer_rle_.n_chars  = rle_start_pos; // reroll back
		return false;
	}

	const U64 n_runs = (this->buffer_rle_.n_chars - rle_start_pos)/this->rleController_->getBitWidth();
	const MetaEntryBase& base_meta = *reinterpret_cast<const MetaEntryBase* const>(&this->buffer_meta_[meta_start_pos]);

	if(n_runs == 1){
		this->buffer_meta_.n_chars = meta_start_pos; // reroll back
		this->buffer_rle_.n_chars  = rle_start_pos; // reroll back
		return false;
	}

	if(base_meta.HWE_P < this->filter.HWE_P){
		this->buffer_meta_.n_chars = meta_start_pos; // reroll back
		this->buffer_rle_.n_chars  = rle_start_pos; // reroll back
		//std::cerr << "HWE_P < " << this->filter.HWE_P << ": " << base_meta.HWE_P << std::endl;
		return false;
	}

	if(base_meta.MAF < this->filter.MAF){
		this->buffer_meta_.n_chars = meta_start_pos; // reroll back
		this->buffer_rle_.n_chars  = rle_start_pos; // reroll back
		//std::cerr << "MAF < " << this->filter.MAF << ": " << base_meta.MAF << std::endl;
		return false;
	}

	if(this->totempole_entry.min_position == 0)
		this->totempole_entry.min_position = line.body->POS + 1;

	this->totempole_entry.max_position = line.body->POS + 1;
	++this->totempole_entry;

	return true;
}

// flush and write
bool ImportWriter::flush(void){
	if(this->buffer_meta_.size() == 0){
		//std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Cannot flush writer with 0 entries..." << std::endl;
		return false;
	}

	this->totempole_entry.byte_offset = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
	this->gzip_controller_.Deflate(this->buffer_meta_, this->buffer_rle_); // Deflate block
	this->streamTomahawk << this->gzip_controller_; // Write tomahawk output
	this->gzip_controller_.Clear(); // Clean up gzip controller

	// Keep track of largest block observed
	if(this->buffer_meta_.size() > this->largest_uncompressed_block_)
		this->largest_uncompressed_block_ = this->buffer_meta_.size();

	this->totempole_entry.uncompressed_size = this->buffer_meta_.size(); // Store uncompressed size
	this->totempole_entry.byte_offset_end = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
	this->streamTotempole << this->totempole_entry; // Write totempole output
	++this->blocksWritten_; // update number of blocks written
	this->variants_written_ += this->totempole_entry.size(); // update number of variants written

	this->reset(); // reset buffers
	return true;
}

void ImportWriter::CheckOutputNames(const std::string& input){
	std::vector<std::string> paths = Helpers::filePathBaseExtension(input);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == Constants::OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &Constants::OUTPUT_SUFFIX[0], Constants::OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}


} /* namespace Tomahawk */
