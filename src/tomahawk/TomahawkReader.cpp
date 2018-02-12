#include "TomahawkReader.h"
#include "genotype_container.h"

namespace Tomahawk {

// Remember to resize buffers to header.getLargestBlockSize()+64 after header is loaded
TomahawkReader::TomahawkReader() :
	filesize_(0),
	offset_end_of_data_(0),
	bit_width_(0),
	dropGenotypes(false),
	showHeader(true),
	index_(nullptr),
	writer(nullptr)
{}

TomahawkReader::~TomahawkReader(){
	this->buffer_.deleteAll();
	this->data_.deleteAll();
	this->outputBuffer_.deleteAll();
	delete this->writer;
	delete this->index_;
}

bool TomahawkReader::open(const std::string input){
	if(input.size() == 0){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "No input filename..." << std::endl;
		return false;
	}

	this->stream_.open(input, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to open file..." << std::endl;
	}
	this->filesize_ = this->stream_.tellg();

	this->stream_.seekg(this->filesize_ - TWK_FOOTER_LENGTH);
	this->stream_ >> this->footer_;
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to open file..." << std::endl;
		return false;
	}

	if(this->footer_.validate() == false){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to validate footer..." << std::endl;
		return false;
	}

	// Seek to start of index
	this->stream_.seekg(this->footer_.offset_end_of_data);
	const U32 l_index_data = (this->filesize_ - TWK_FOOTER_LENGTH) - this->stream_.tellg();
	buffer_type index_buffer(l_index_data + 1024);
	this->stream_.read(index_buffer.data(), l_index_data);
	index_buffer.n_chars = l_index_data;
	this->index_ = new index_type(index_buffer.data(), index_buffer.size());
	index_buffer.deleteAll();

	// Resize buffers to accomodate the largest possible block
	// without ever resizing
	// this is for performance reasons
	this->buffer_.resize(this->getFooter().getLargestUncompressedBlock() + 64);
	this->data_.resize(this->getFooter().getLargestUncompressedBlock() + 64);
	this->outputBuffer_.resize(this->getFooter().getLargestUncompressedBlock() + 64);

	// Seek to beginning
	this->stream_.seekg(0);
	if(!this->header_.open(this->stream_)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to load header data..." << std::endl;
		return false;
	}

	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Stream is bad..." << std::endl;
		return false;
	}

	if(this->header_.validate() == false){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to validate header..." << std::endl;
		return false;
	}

	this->offset_end_of_data_ = this->footer_.offset_end_of_data;
	this->DetermineBitWidth();

	return true;
}

bool TomahawkReader::getBlocks(void){
	U64 buffer_size = 0;
	for(U32 i = 0; i < this->index_->getContainer().size(); ++i){
		buffer_size += this->index_->getContainer()[i].uncompressed_size;
	}

	if(buffer_size == 0){
		std::cerr << "impsosible" << std::endl;
		return false;
	}

	//if(!SILENT)
	//	std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << this->totempole_.getHeader().getNumberBlocks() << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

	this->data_.resize(buffer_size + 1000);

	for(U32 i = 0; i < this->index_->getContainer().size(); ++i){
		if(!this->getBlock(i))
			return false;
	}

	return true;
}

bool TomahawkReader::getBlocks(std::vector<U32>& blocks){
	U64 buffer_size = 0;
	for(U32 i = 0; i < blocks.size(); ++i){
		std::cerr << i << '/' << blocks.size() << '\t' << this->index_->getContainer()[i].uncompressed_size << std::endl;
		buffer_size += this->index_->getContainer()[i].uncompressed_size;
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
			buffer_size += this->index_->getContainer()[j].uncompressed_size;
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

	this->stream_.seekg(this->index_->getContainer()[blockID].byte_offset);
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed search..." << std::endl;
		return false;
	}

	const U32 readLength = this->index_->getContainer()[blockID].byte_offset_end - this->index_->getContainer()[blockID].byte_offset;

	if(readLength > this->buffer_.capacity()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Overflowing capacity: " << readLength << '/' << this->buffer_.capacity() << std::endl;
		exit(1);
	}

	this->blockDataOffsets_.push_back(DataOffsetPair(&this->data_[this->data_.size()], this->index_->getContainer()[blockID].uncompressed_size, this->index_->getContainer()[blockID]));
	if(!this->stream_.read(this->buffer_.data(), readLength)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed read: " << this->stream_.good() << '\t' << this->stream_.fail() << '/' << this->stream_.eof() << std::endl;
		//std::cerr << this->stream_.gcount() << '/' << readLength << std::endl;
		return false;
	}
	this->buffer_.n_chars = readLength;

	if(!this->tgzf_controller_.Inflate(this->buffer_, this->data_)){
		std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to inflate data..." << std::endl;
		return false;
	}

	return true;
}

void TomahawkReader::DetermineBitWidth(void){
	if(this->header_.magic_.getNumberSamples() <= Constants::UPPER_LIMIT_SAMPLES_8B - 1){
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " < " << Constants::UPPER_LIMIT_SAMPLES_8B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 8-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(BYTE);
	} else if(this->header_.magic_.getNumberSamples() <= Constants::UPPER_LIMIT_SAMPLES_16B - 1){
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " < " << Constants::UPPER_LIMIT_SAMPLES_16B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 16-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(U16);
	} else if(this->header_.magic_.getNumberSamples() <= Constants::UPPER_LIMIT_SAMPLES_32B - 1){
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " < " << Constants::UPPER_LIMIT_SAMPLES_32B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 32-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(U32);
	} else {
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << Constants::UPPER_LIMIT_SAMPLES_32B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " < " << Constants::UPPER_LIMIT_SAMPLES_64B << "..." << std::endl;
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
		std::cout << this->header_.getLiterals() + "\n##tomahawk_viewCommand=" + Helpers::program_string() << std::endl;


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
		std::cerr << Helpers::timestamp("LOG", "TGZF") << "Inflating " << this->index_->getContainer().size() << " blocks..." << std::endl;

	// Output header
	if(this->showHeader){
		std::cout << this->header_.getLiterals() << std::endl;
		std::cout << "##INFO=<ID=HWE_P,Number=1,Type=Float,Description=\"Hardy-Weinberg P-value (Fisher's exact test)\">" << std::endl;
		std::cout << "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">" << std::endl;
		std::cout << "##tomahawk_viewCommand=" + Helpers::program_string() << std::endl;
		std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";

		std::cout << this->header_.getSample(0);
		for(U32 i = 1; i < this->header_.magic_.getNumberSamples(); ++i)
			std::cout << '\t' << this->header_.getSample(i);
		std::cout.put('\n');

	}

	for(U32 i = 0; i < this->index_->getContainer().size(); ++i){
		(*this.*func__)(i);
	}

	return true;
}


} /* namespace Tomahawk */
