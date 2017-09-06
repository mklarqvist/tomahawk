#include "TomahawkReader.h"

namespace Tomahawk {

// Remember to resize buffers to header.getLargestBlockSize()+64 after header is loaded
TomahawkReader::TomahawkReader() :
	samples(0),
	version(0),
	filesize_(0),
	bit_width_(0),
	dropGenotypes(false),
	showHeader(true),
	writer(nullptr)
{}

TomahawkReader::~TomahawkReader(){
	this->buffer_.deleteAll();
	this->data_.deleteAll();
	this->outputBuffer_.deleteAll();
	delete writer;
}

bool TomahawkReader::Open(const std::string input){
	const std::string index = input + '.' + Tomahawk::Constants::OUTPUT_INDEX_SUFFIX;

	// Parse Totempole
	if(!this->totempole_.Open(index)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TOTEMPOLE") << "Failed build!" << std::endl;
		return false;
	}

	// Resize buffers to accomodate the largest possible block
	// without ever resizing
	// this is for performance reasons
	this->buffer_.resize(this->totempole_.getLargestBlockSize() + 64);
	this->data_.resize(this->totempole_.getLargestBlockSize() + 64);
	this->outputBuffer_.resize(this->totempole_.getLargestBlockSize() + 64);

	if(input.size() == 0){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "No input filename..." << std::endl;
		return false;
	}

	this->stream_.open(input, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Could not open " << input << "..." << std::endl;
		return false;
	}
	this->filesize_ = this->stream_.tellg();
	this->stream_.seekg(0);

	// Validate MAGIC and header
	if(!this->ValidateHeader(this->stream_)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to validate header..." << std::endl;
		return false;
	}

	return(this->Validate());
}


inline bool TomahawkReader::ValidateHeader(std::ifstream& in) const{
	char MAGIC[Constants::WRITE_HEADER_MAGIC_LENGTH];
	in.read(MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH);

	if(strncmp(&MAGIC[0], &Constants::WRITE_HEADER_MAGIC[0], Constants::WRITE_HEADER_MAGIC_LENGTH) == 0)
		return true;

	return false;
}

bool TomahawkReader::getBlocks(void){
	U64 buffer_size = 0;
	for(U32 i = 0; i < this->totempole_.header.blocks; ++i){
		buffer_size += this->totempole_[i].uncompressed_size;
	}

	if(buffer_size == 0){
		std::cerr << "impsosible" << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << this->totempole_.getHeader().blocks << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

	this->data_.resize(buffer_size + 1000);

	for(U32 i = 0; i < this->totempole_.getBlocks(); ++i)
		if(!this->getBlock(i))
			return false;

	return true;
}

bool TomahawkReader::getBlocks(std::vector<U32>& blocks){
	U64 buffer_size = 0;
	for(U32 i = 0; i < blocks.size(); ++i){
		std::cerr << i << '/' << blocks.size() << '\t' << this->totempole_[i].uncompressed_size << std::endl;
		buffer_size += this->totempole_[i].uncompressed_size;
	}

	if(buffer_size == 0)
		return false;

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << blocks.size() << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

	this->data_.resize(buffer_size + 1000);

	for(U32 i = 0; i < blocks.size(); ++i)
		if(!this->getBlock(blocks[i]))
			return false;

	return true;
}

bool TomahawkReader::getBlocks(std::vector< std::pair<U32, U32> >& blocks){
	U64 buffer_size = 0;
	for(U32 i = 0; i < blocks.size(); ++i){
		for(U32 j = blocks[i].first; j < blocks[i].second; ++j){
			buffer_size += this->totempole_[j].uncompressed_size;
		}
	}

	if(buffer_size == 0){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Uncompressed size is 0..." << std::endl;
		return false;
	}

	U64 totalBlocks = 0;
	for(U32 i = 0; i < blocks.size(); ++i)
		totalBlocks += blocks[i].second - blocks[i].first;

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << totalBlocks << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

	this->data_.resize(buffer_size + 1000);

	for(U32 i = 0; i < blocks.size(); ++i){
		for(U32 j = blocks[i].first; j < blocks[i].second; ++j){
			if(!this->getBlock(j)){
				std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Could not get block: " << j << '/' << blocks[i].second << std::endl;
				return false;
			}
		}
	}
	return true;
}


bool TomahawkReader::getBlock(const U32 blockID){
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Bad file stream for block " << blockID << "..." << std::endl;
		return false;
	}

	this->stream_.seekg(this->totempole_[blockID].byte_offset);
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed search..." << std::endl;
		return false;
	}

	const U32 readLength = this->totempole_[blockID].byte_offset_end - this->totempole_[blockID].byte_offset;

	if(readLength > this->buffer_.capacity()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Overflowing capacity: " << readLength << '/' << this->buffer_.capacity() << std::endl;
		exit(1);
	}

	this->blockDataOffsets_.push_back(DataOffsetPair(&this->data_.data[this->data_.pointer], this->totempole_[blockID]));
	if(!this->stream_.read(&this->buffer_.data[0], readLength)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed read: " << this->stream_.good() << '\t' << this->stream_.fail() << '/' << this->stream_.eof() << std::endl;
		//std::cerr << this->stream_.gcount() << '/' << readLength << std::endl;
		return false;
	}
	this->buffer_.pointer = readLength;

	if(!this->tgzf_controller_.Inflate(this->buffer_, this->data_)){
		std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to inflate data..." << std::endl;
		return false;
	}

	return true;
}

bool TomahawkReader::Validate(void){
	char temp_buffer[sizeof(float)+sizeof(U64)];
	this->stream_.read(&temp_buffer[0], sizeof(float)+sizeof(U64));

	const float* version = reinterpret_cast<const float*>(&temp_buffer[0]);
	const U64* samples   = reinterpret_cast<const U64*>(&temp_buffer[sizeof(float)]);

	this->version = *version;
	this->samples = *samples;

	if(this->version != this->totempole_.getHeader().version){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "File discordance: versions do not match..." << std::endl;
		return false;
	}

	if(this->samples != this->totempole_.getHeader().samples){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "File discordance:number of samples do not match" << std::endl;
		return false;
	}

	// Determine what bit-width functions to use
	this->DetermineBitWidth();

	return true;
}

void TomahawkReader::DetermineBitWidth(void){
	if(this->samples <= Constants::UPPER_LIMIT_SAMPLES_8B - 1){
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_8B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 8-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(BYTE);
	} else if(this->samples <= Constants::UPPER_LIMIT_SAMPLES_16B - 1){
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_16B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 16-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(U16);
	} else if(this->samples <= Constants::UPPER_LIMIT_SAMPLES_32B - 1){
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_32B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 32-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(U32);
	} else {
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_32B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_64B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 64-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(U64);
	}
}

bool TomahawkReader::outputBlocks(std::vector<U32>& blocks){
	typedef bool (Tomahawk::TomahawkReader::*blockFunction)(const U32 blockID);
	blockFunction func__ = nullptr;

	switch(this->bit_width_){
	case 1: func__ = &TomahawkReader::outputBlock<BYTE>; break;
	case 2: func__ = &TomahawkReader::outputBlock<U16>;  break;
	case 4: func__ = &TomahawkReader::outputBlock<U32>;  break;
	case 8: func__ = &TomahawkReader::outputBlock<U64>;  break;
	default: exit(1); break;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","TGZF") << "Inflating " << blocks.size() << " blocks..." << std::endl;

	// Output header
	if(this->showHeader)
		std::cout << this->totempole_.literals + "\n##tomahawk_viewCommand=" + Helpers::program_string(true) << std::endl;


	for(U32 i = 0; i < blocks.size(); ++i)
		(*this.*func__)(blocks[i]);

	return true;
}

bool TomahawkReader::outputBlocks(){
	typedef bool (Tomahawk::TomahawkReader::*blockFunction)(const U32 blockID);
	blockFunction func__ = nullptr;

	switch(this->bit_width_){
	case 1: func__ = &TomahawkReader::outputBlock<BYTE>; break;
	case 2: func__ = &TomahawkReader::outputBlock<U16>;  break;
	case 4: func__ = &TomahawkReader::outputBlock<U32>;  break;
	case 8: func__ = &TomahawkReader::outputBlock<U64>;  break;
	default:
		std::cerr << Helpers::timestamp("ERROR") << "Word sizing could not be determined!" << std::endl;
		exit(1);
		break;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "TGZF") << "Inflating " << this->totempole_.getBlocks() << " blocks..." << std::endl;

	// Output header
	if(this->showHeader){
		std::cout << this->totempole_.literals << std::endl;
		std::cout << "##INFO=<ID=HWE_P,Number=1,Type=Float,Description=\"Hardy-Weinberg P-value (Fisher's exact test)\">" << std::endl;
		std::cout << "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">" << std::endl;
		std::cout << "##tomahawk_viewCommand=" + Helpers::program_string(true) << std::endl;
		std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
		for(U32 i = 0; i < this->totempole_.header.samples - 1; ++i)
			std::cout << this->totempole_.samples[i] << '\t';
		std::cout << this->totempole_.samples[this->totempole_.header.samples - 1] << std::endl;
	}

	for(U32 i = 0; i < this->totempole_.getBlocks(); ++i)
		(*this.*func__)(i);

	return true;
}


} /* namespace Tomahawk */
