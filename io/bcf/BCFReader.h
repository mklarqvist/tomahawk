
#ifndef BCFREADER_H_
#define BCFREADER_H_

#include "../BasicBuffer.h"
#include "../BGZFController.h"

namespace Tomahawk {
namespace IO {

#pragma pack(1)
struct BCFEntryBody{
	BCFEntryBody(); // disallow ctor and dtor
	~BCFEntryBody();

	U32 l_shared;
	U32 l_indiv;
	S32 CHROM;
	S32 POS;
	S32 rlen;
	float QUAL;
	U32 n_allele: 16, n_info: 16;
	U32 n_fmt: 24, n_sample: 8;
};

class BCFReader{
	typedef BCFReader self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::BGZFController bgzf_controller_type;
	typedef BGZFHeader bgzf_type;

public:
	BCFReader(){}
	BCFReader(const std::string file);
	~BCFReader(){}

	bool nextBlock(void){
		// Stream died
		if(!this->stream.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Stream died!" << std::endl;
			return false;
		}

		// EOF
		if(this->stream.tellg() == this->filesize){
			//std::cerr << "eof" << std::endl;
			return false;
		}

		std::cerr << this->stream.tellg() << std::endl;
		buffer.resize(sizeof(bgzf_type));
		this->stream.read(&buffer.data[0],  Constants::TGZF_BLOCK_HEADER_LENGTH);
		const bgzf_type* h = reinterpret_cast<const bgzf_type*>(&buffer.data[0]);
		buffer.pointer = Constants::TGZF_BLOCK_HEADER_LENGTH;
		if(!h->Validate()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Failed to validate!" << std::endl;
			std::cerr << *h << std::endl;
			return false;
		}

		buffer.resize(h->BSIZE); // make sure all data will fit

		// Recast because if buffer is resized then the pointer address is incorrect
		// resulting in segfault
		h = reinterpret_cast<const bgzf_type*>(&buffer.data[0]);

		this->stream.read(&buffer.data[Constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - Constants::TGZF_BLOCK_HEADER_LENGTH);
		if(!this->stream.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Truncated file..." << std::endl;
			return false;
		}

		buffer.pointer = h->BSIZE;
		const U32 uncompressed_size = *reinterpret_cast<const U32*>(&buffer[buffer.pointer -  sizeof(U32)]);
		output_buffer.resize(uncompressed_size);
		this->output_buffer.reset();

		if(!this->bgzf_controller.Inflate(buffer, output_buffer)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Failed inflate!" << std::endl;
			return false;
		}

		if(this->output_buffer.size() == 0){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Empty data!" << std::endl;
			return false;
		}

		// Reset buffer
		this->buffer.reset();

		return true;
	}

	bool parseHeader(void){
		if(this->output_buffer.size() == 0){
			std::cerr << "no buffer" << std::endl;
			return false;
		}

		std::cerr << "in parse header" << std::endl;

		std::cerr << std::string(&this->output_buffer[0], 5) << std::endl;
		const U32& l_text = *reinterpret_cast<const U32* const>(&this->output_buffer[5]);
		std::cerr << "length: " << l_text << std::endl;
		return true;
	}

	bool open(const std::string input){
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

private:
	std::ifstream stream;
	U64 filesize;
	buffer_type buffer;
	buffer_type output_buffer;
	bgzf_controller_type bgzf_controller;
};


}
}



#endif /* BCFREADER_H_ */
