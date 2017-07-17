#include "TomahawkReader.h"

namespace Tomahawk {

bool TomahawkReader::Open(const std::string input){
	if(input.size() == 0){
		std::cerr << "no filename" << std::endl;
		return false;
	}

	this->stream_.open(input, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream_.good()){
		std::cerr << "could not open " << input << "..." << std::endl;
		return false;
	}
	this->filesize_ = this->stream_.tellg();
	this->stream_.seekg(0);

	return true;
}

bool TomahawkReader::ValidateHeader(void){
	if(!this->stream_.good()){
		std::cerr << "invalid filehandle..." << std::endl;
		return false;
	}

	// Read some data
	this->stream_.read(this->buffer_.data, Constants::TOMAHAWK_HEADER_LENGTH);
	this->buffer_.pointer += Constants::TOMAHAWK_HEADER_LENGTH;

	return this->Validate();
}

U64 TomahawkReader::GetUncompressedSizes(std::vector< std::pair<U32, U32> >& blocks){
	if(blocks.size() == 0)
		return 0;

	char temp[sizeof(U32)];
	U32* BSIZE = reinterpret_cast<U32*>(&temp[0]);
	U64 totalSize = 0;

	U32 number_of_blocks = blocks.size();
	bool adjusted = false;

	if(blocks[blocks.size()-1].second == this->totempole_.getBlocks()){
		//std::cerr << "Last one is final" << std::endl;
		this->stream_.seekg(this->filesize_ - Constants::eof_length*sizeof(U64) - sizeof(U32));
		if(!this->stream_.good()){
			std::cerr << "Illegal seek" << std::endl;
			return(0);
		}
		this->stream_.read(&temp[0], sizeof(U32));
		totalSize += *BSIZE;
		adjusted = true;
		--blocks[blocks.size()-1].second;
	}
	// Check for duplicates

	for(U32 i = 0; i < number_of_blocks; ++i){
		for(U32 j = blocks[i].first; j < blocks[i].second; ++j){
			this->stream_.seekg(this->totempole_[j+1].byte_offset - sizeof(U32));
			if(!this->stream_.good()){
				std::cerr << "Illegal seek" << std::endl;
				return(0);
			}
			this->stream_.read(&temp[0], sizeof(U32));
			//std::cerr << *BSIZE << std::endl;
			totalSize += *BSIZE;
		}
	}

	// Restore if last one was last
	if(adjusted)
		++blocks[blocks.size()-1].second;

	return totalSize;
}

U64 TomahawkReader::GetUncompressedSizes(std::vector<U32>& blocks){
	if(blocks.size() == 0)
		return 0;

	char temp[sizeof(U32)];
	U32* BSIZE = reinterpret_cast<U32*>(&temp[0]);
	U64 totalSize = 0;

	U32 number_of_blocks = blocks.size();

	std::sort(blocks.begin(), blocks.end());
	if(blocks[blocks.size()-1] == this->totempole_.getBlocks() -1){
		//std::cerr << "Last one is final" << std::endl;
		this->stream_.seekg(this->filesize_ - Constants::eof_length*sizeof(U64) - sizeof(U32));
		if(!this->stream_.good()){
			std::cerr << "Illegal seek" << std::endl;
			return(0);
		}
		this->stream_.read(&temp[0], sizeof(U32));
		totalSize += *BSIZE;
		--number_of_blocks;
	}
	// Check for duplicates

	for(U32 i = 0; i < number_of_blocks; ++i){
		this->stream_.seekg(this->totempole_[blocks[i]+1].byte_offset - sizeof(U32));
		if(!this->stream_.good()){
			std::cerr << "Illegal seek" << std::endl;
			return(0);
		}
		this->stream_.read(&temp[0], sizeof(U32));
		//std::cerr << *BSIZE << std::endl;
		totalSize += *BSIZE;
	}

	return totalSize;
}

U64 TomahawkReader::GetUncompressedSizes(void){
	char temp[sizeof(U32)];
	U32* BSIZE = reinterpret_cast<U32*>(&temp[0]);
	U64 totalSize = 0;

	for(U32 i = 0; i < this->totempole_.getBlocks()-1; ++i){
		this->stream_.seekg(this->totempole_[i+1].byte_offset - sizeof(U32));
		this->stream_.read(&temp[0], sizeof(U32));
		if(!this->stream_.good()){
			std::cerr << "Illegal seek" << std::endl;
			return(0);
		}
		totalSize += *BSIZE;
	}

	this->stream_.seekg(this->filesize_ - Constants::eof_length*sizeof(U64) - sizeof(U32));
	if(!this->stream_.good()){
		std::cerr << "Illegal seek" << std::endl;
		return(0);
	}
	this->stream_.read(&temp[0], sizeof(U32));
	totalSize += *BSIZE;

	//std::cerr << "Final: " << *BSIZE << std::endl;

	return totalSize;
}

bool TomahawkReader::getBlocks(void){
	U64 buffer_size = this->GetUncompressedSizes();
	if(buffer_size == 0)
		return false;

	std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << this->totempole_.getHeader().blocks << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;
	this->data_.resize(buffer_size + 1000);

	for(U32 i = 0; i < this->totempole_.getBlocks(); ++i)
		if(!this->getBlock(i))
			return false;

	//std::cerr << "Data: " << this->data_.size() << " expected " << buffer_size << std::endl;

	return true;
}

bool TomahawkReader::getBlocks(std::vector<U32>& blocks){
	U64 buffer_size = this->GetUncompressedSizes(blocks);
	if(buffer_size == 0)
		return false;

	std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << blocks.size() << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;
	this->data_.resize(buffer_size + 1000);

	for(U32 i = 0; i < blocks.size(); ++i)
		if(!this->getBlock(blocks[i]))
			return false;

	//std::cerr << "Data: " << this->data_.size() << " expected " << buffer_size << std::endl;

	return true;
}

bool TomahawkReader::getBlocks(std::vector< std::pair<U32, U32> >& blocks){
	U64 buffer_size = this->GetUncompressedSizes(blocks);
	if(buffer_size == 0){
		std::cerr << "Buffer size is 0" << std::endl;
		return false;
	}

	U64 totalBlocks = 0;
	for(U32 i = 0; i < blocks.size(); ++i)
		totalBlocks += blocks[i].second - blocks[i].first;

	std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << totalBlocks << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;
	this->data_.resize(buffer_size + 1000);

	for(U32 i = 0; i < blocks.size(); ++i){
		for(U32 j = blocks[i].first; j < blocks[i].second; ++j){
			if(!this->getBlock(j)){
				std::cerr << "could not get block: " << j << '/' << blocks[i].second << std::endl;
				return false;
			}
		}
	}
	return true;
}


bool TomahawkReader::getBlock(const U32 blockID){
	if(!this->stream_.good()){
		std::cerr << "stream bad " << blockID << std::endl;
		return false;
	}

	//std::cerr << "getblock " << blockID  << " seek to " << this->totempole_[blockID].byte_offset << std::endl;
	this->stream_.seekg(this->totempole_[blockID].byte_offset);
	if(!this->stream_.good()){
		std::cerr << "Failed search" << std::endl;
		return false;
	}

	U32 readLength = 0;
	if(blockID != this->totempole_.getBlocks() - 1)
		readLength = this->totempole_[blockID + 1].byte_offset - this->totempole_[blockID].byte_offset;
	else
		readLength = this->filesize_ - Constants::eof_length*sizeof(U64) - this->totempole_[this->totempole_.getBlocks()-1].byte_offset;

	if(readLength > this->buffer_.capacity()){
		std::cerr << "impossible: " << readLength << '/' << this->buffer_.capacity() << std::endl;
		exit(1);
	}

	//this->manager_.Add(&this->data_.data[this->data_.pointer], this->totempole_[blockID]);
	this->blockDataOffsets_.push_back(DataOffsetPair(&this->data_.data[this->data_.pointer], this->totempole_[blockID]));
	if(!this->stream_.read(&this->buffer_.data[0], readLength)){
		std::cerr << "Failed read: " << this->stream_.good() << '\t' << this->stream_.fail() << '/' << this->stream_.eof() << std::endl;
		std::cerr << this->stream_.gcount() << '/' << readLength << std::endl;
		return false;
	}
	this->buffer_.pointer = readLength;

	if(!this->gzip_controller_.Inflate(this->buffer_, this->data_)){
		std::cerr << "failed" << std::endl;
		return false;
	}

	//std::cerr << "Done: " << blockID << std::endl;

	return true;
}

bool TomahawkReader::Validate(void){
	if(this->buffer_.size() < Constants::WRITE_HEADER_MAGIC_LENGTH + sizeof(float) + sizeof(U64)){
		std::cerr << "failed to validate: corrupt" << std::endl;
		return false;
	}

	if(strncmp(this->buffer_.data, Constants::WRITE_HEADER_MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH) != 0){
		std::cerr << "illegal tomahawk file" << std::endl;
		return false;
	}

	const float* version = reinterpret_cast<const float*>(&this->buffer_[Constants::WRITE_HEADER_MAGIC_LENGTH]);
	const U64* samples = reinterpret_cast<const U64*>(&this->buffer_[Constants::WRITE_HEADER_MAGIC_LENGTH + sizeof(float)]);

	this->version = *version;
	this->samples = *samples;
	this->progress.SetSamples(this->samples); // push to progressbar

	if(this->version != this->totempole_.getHeader().version){
		std::cerr << "versions do not match" << std::endl;
		return false;
	}

	if(this->samples != this->totempole_.getHeader().samples){
		std::cerr << "number of samples do not match" << std::endl;
		return false;
	}

	// Determine what bit-width functions to use
	this->DetermineBitWidth();

	return true;
}

void TomahawkReader::DetermineBitWidth(void){
	if(this->samples <= Constants::UPPER_LIMIT_SAMPLES_8B - 1){
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_8B << "..." << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 8-bit width..." << std::endl;
		this->bit_width_ = sizeof(BYTE);
	} else if(this->samples <= Constants::UPPER_LIMIT_SAMPLES_16B - 1){
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_16B << "..." << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 16-bit width..." << std::endl;
		this->bit_width_ = sizeof(U16);
	} else if(this->samples <= Constants::UPPER_LIMIT_SAMPLES_32B - 1){
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_32B << "..." << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 32-bit width..." << std::endl;
		this->bit_width_ = sizeof(U32);
	} else {
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_32B << "... Skip" << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_64B << "..." << std::endl;
		std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 64-bit width..." << std::endl;
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

	std::cerr << Helpers::timestamp("LOG","BGZF") << "Inflating " << blocks.size() << " blocks..." << std::endl;
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
	default: exit(1); break;
	}

	std::cerr << Helpers::timestamp("LOG", "BGZF") << "Inflating " << this->totempole_.getBlocks() << " blocks..." << std::endl;
	for(U32 i = 0; i < this->totempole_.getBlocks(); ++i)
		(*this.*func__)(i);

	return true;
}

template <typename K, typename V>
bool comparePairs(const std::pair<K,V>& a, const std::pair<K,V>& b){ return a.first < b.first; }

bool TomahawkReader::__CalculateWrapper(){
	if(this->bit_width_ == 1) 	   return(this->__Calculate<BYTE>());
	else if(this->bit_width_ == 2) return(this->__Calculate<U16>());
	else if(this->bit_width_ == 4) return(this->__Calculate<U32>());
	else if(this->bit_width_ == 8) return(this->__Calculate<U64>());
	else {
		std::cerr << "impossible bit width" << std::endl;
		exit(1);
	}

	//std::cerr << "Manager have: " << this->manager_.size() << std::endl;
	//return(this->manager_.AllVersusAll());
	return false;
}

bool TomahawkReader::Calculate(std::vector< std::pair<U32,U32> >& blocks){
	// Todo!
	std::sort(blocks.begin(), blocks.end(), comparePairs<U32, U32>);
	if(!this->getBlocks(blocks)){
		std::cerr << "Failed to get blocks" << std::endl;
		return false;
	}

	return(this->__CalculateWrapper());
}

bool TomahawkReader::Calculate(std::vector<U32>& blocks){
	if(!this->getBlocks(blocks)){
		std::cerr << "Failed" << std::endl;
		return false;
	}

	std::cerr << Helpers::timestamp("LOG") << "Inflated " << blocks.size() << " blocks..." << std::endl;
	return(this->__CalculateWrapper());
}

bool TomahawkReader::Calculate(){
	if(!this->balancer.Build(this->totempole_.getBlocks(), this->threads)){
		std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Failed to split into blocks..." << std::endl;
		return false;
	}

	return(this->Calculate(this->balancer.getLoad()));
}

bool TomahawkReader::SelectWriterOutputType(const IO::TomahawkCalculationWriterInterace::type& writer_type){
	if(this->writer != nullptr)
		return false;

	if(writer_type == IO::TomahawkCalculationWriterInterace::type::cout)
		this->writer = new IO::TomahawkCalculationWriterStandardOut;
	else
		this->writer = new IO::TomahawkCalculationWriterFile;

	return true;
}

} /* namespace Tomahawk */
