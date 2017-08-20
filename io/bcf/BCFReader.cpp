#include "../vcf/VCFHeader.h"
#include "BCFReader.h"

#ifndef IO_BCF_BCFREADER_CPP_
#define IO_BCF_BCFREADER_CPP_

namespace Tomahawk{
namespace BCF{

BCFReader::BCFReader() : filesize(0), current_pointer(0){}
BCFReader::~BCFReader(){}


bool BCFReader::nextBlock(void){
	// Stream died
	if(!this->stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Stream died!" << std::endl;
		return false;
	}

	// EOF
	if(this->stream.tellg() == this->filesize)
		return false;

	buffer.resize(sizeof(bgzf_type));
	this->stream.read(&buffer.data[0], IO::Constants::BGZF_BLOCK_HEADER_LENGTH);
	const bgzf_type* h = reinterpret_cast<const bgzf_type*>(&buffer.data[0]);
	buffer.pointer = IO::Constants::BGZF_BLOCK_HEADER_LENGTH;
	if(!h->Validate()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Failed to validate!" << std::endl;
		std::cerr << *h << std::endl;
		return false;
	}

	buffer.resize(h->BSIZE); // make sure all data will fit

	// Recast because if buffer is resized then the pointer address is incorrect
	// resulting in segfault
	h = reinterpret_cast<const bgzf_type*>(&buffer.data[0]);

	this->stream.read(&buffer.data[IO::Constants::BGZF_BLOCK_HEADER_LENGTH], (h->BSIZE + 1) - IO::Constants::BGZF_BLOCK_HEADER_LENGTH);
	if(!this->stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Truncated file..." << std::endl;
		return false;
	}

	buffer.pointer = h->BSIZE + 1;
	const U32 uncompressed_size = *reinterpret_cast<const U32*>(&buffer[buffer.pointer -  sizeof(U32)]);
	output_buffer.resize(uncompressed_size);
	this->output_buffer.reset();

	if(!this->bgzf_controller.Inflate(buffer, output_buffer)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Failed inflate!" << std::endl;
		return false;
	}

	// BGZF EOF marker
	if(this->output_buffer.size() == 0)
		return false;

	// Reset buffer
	this->buffer.reset();
	this->current_pointer = 0;

	return true;
}

bool BCFReader::nextVariant(BCFEntry& entry){
	if(this->current_pointer + 8 > this->output_buffer.size()){
		const U32 partial = this->output_buffer.size() - this->current_pointer;
		entry.add(&this->output_buffer[this->current_pointer], this->output_buffer.size() - this->current_pointer);
		if(!this->nextBlock()){
			std::cerr << "failed to get next block" << std::endl;
			return false;
		}

		entry.add(&this->output_buffer[0], 8 - partial);
		this->current_pointer = 8 - partial;
	} else {
		entry.add(&this->output_buffer[this->current_pointer], 8);
		this->current_pointer += 8;
	}

	U64 remainder = entry.sizeBody();
	while(remainder > 0){
		if(this->current_pointer + remainder > this->output_buffer.size()){
			entry.add(&this->output_buffer[this->current_pointer], this->output_buffer.size() - this->current_pointer);
			remainder -= this->output_buffer.size() - this->current_pointer;
			if(!this->nextBlock()){
				std::cerr << "failed to get next block" << std::endl;
				return false;
			}
		} else {
			entry.add(&this->output_buffer[this->current_pointer], remainder);
			this->current_pointer += remainder;
			remainder = 0;
			break;
		}
	}

	// Temp
	//std::cerr << this->header.getContig(entry.body->CHROM).name << ":" << entry.body->POS << std::endl;
	entry.parse();

	return true;
}

bool BCFReader::parseHeader(void){
	if(this->output_buffer.size() == 0){
		std::cerr << "no buffer" << std::endl;
		return false;
	}

	if(strncmp(&this->output_buffer.data[0], "BCF\2\2", 5) != 0){ // weird: should be BCF/2/1
		std::cerr << (int)this->output_buffer[3] << '\t' << (int)this->output_buffer[4] << std::endl;
		std::cerr << "failed to validate" << std::endl;
		return false;
	}

	const U32 l_text = *reinterpret_cast<const U32* const>(&this->output_buffer[5]) + 4;
	this->header_buffer.resize(l_text);

	if(l_text - 5 < this->output_buffer.size()){
		this->header_buffer.Add(&this->output_buffer[5], l_text);
		return true;
	}

	U32 head_read = this->output_buffer.size() - 5;
	this->header_buffer.Add(&this->output_buffer[5], this->output_buffer.size() - 5);

	//U32 p = 0;
	while(this->nextBlock()){
		if(head_read + this->output_buffer.size() >= l_text){
			std::cerr << "remainder: " << l_text - head_read << " and data: " << this->output_buffer.size() << std::endl;
			this->header_buffer.Add(&this->output_buffer[0], l_text - head_read);
			this->current_pointer = l_text - head_read;
			break;
		}
		head_read += this->output_buffer.size();
		this->header_buffer.Add(&this->output_buffer[0], this->output_buffer.size());
	}


	if(!this->header.parse(&this->header_buffer[0], this->header_buffer.size())){
		std::cerr << "failed to parse header" << std::endl;
		return false;
	}

	return true;
}

bool BCFReader::open(const std::string input){
	this->stream.open(input, std::ios::binary | std::ios::in | std::ios::ate);
	if(!this->stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Failed to open file: " << input << std::endl;
		return false;
	}

	this->filesize = this->stream.tellg();
	this->stream.seekg(0);

	if(!this->stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Bad stream!" << std::endl;
		return false;
	}

	if(!this->nextBlock()){
		std::cerr << "failed ot get first block" << std::endl;
		return false;
	}

	if(!this->parseHeader()){
		std::cerr << "failed to parse bcf header" << std::endl;
		return false;
	}

	return true;
}

}
}



#endif /* IO_BCF_BCFREADER_CPP_ */
