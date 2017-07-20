#include "TomahawkReader.h"

namespace Tomahawk {

TomahawkReader::TomahawkReader(const TotempoleReader& header) :
	samples(0),
	version(0),
	filesize_(0),
	bit_width_(0),
	threads(std::thread::hardware_concurrency() > 0 ? std::thread::hardware_concurrency() : 1),
	writer(nullptr),
	buffer_(header.getLargestBlockSize()+64),
	data_(this->buffer_),
	outputBuffer_(header.getLargestBlockSize()+64),
	gzip_controller_(),
	totempole_(header),
	balancer()
{}

TomahawkReader::~TomahawkReader(){
	this->buffer_.deleteAll();
	this->data_.deleteAll();
	this->outputBuffer_.deleteAll();
	delete writer;
}

bool TomahawkReader::SetR2Threshold(const double min, const double max){
		if(min < 0 || max > 1){
			std::cerr << Helpers::timestamp("ERROR") << "Invalid R-squared values: " << min << '-' << max << std::endl;
			return false;
		}

		if(min > max){
			std::cerr << Helpers::timestamp("ERROR") << "Invalid R-squared values: " << min << '-' << max << std::endl;
			return false;
		}

		this->parameters.R2_min = min - Constants::ALLOWED_ROUNDING_ERROR;
		this->parameters.R2_max = max + Constants::ALLOWED_ROUNDING_ERROR;

		return true;
	}

bool TomahawkReader::SetMinimumAlleles(const U64 min){
	this->parameters.minimum_alleles = min;
	return true;
}

bool TomahawkReader::SetThreads(const S32 threads){
	if(threads <= 0){
		std::cerr << Helpers::timestamp("ERROR") << "Invalid number of threads: " << threads << std::endl;
		return false;
	}
	this->threads = threads;
	return true;
}

void TomahawkReader::SetPhased(const bool phased){
	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG") << "Forcing phasing of all variants: " << (phased ? "Phased..." : "Unphased...") << std::endl;

	if(phased)
		this->parameters.force = parameter_type::phasedFunction;
	else
		this->parameters.force = parameter_type::unphasedFunction;
}

bool TomahawkReader::SetPThreshold(const double P){
		if(P > 1 || P < 0){
			std::cerr << Helpers::timestamp("ERROR") << "Invalid P-value threshold: " << P << std::endl;
			return false;
		}
		this->parameters.P_threshold = P;
		return true;
	}

bool TomahawkReader::OpenWriter(void){
	if(this->writer == nullptr)
		this->SelectWriterOutputType(IO::GenericWriterInterace::type::cout);

	return(this->writer->open());
}


bool TomahawkReader::OpenWriter(const std::string destination){
	if(this->writer == nullptr){
		if(destination == "-"){
			this->SelectWriterOutputType(IO::GenericWriterInterace::type::cout);
			return(this->writer->open());
		}
		this->SelectWriterOutputType(IO::GenericWriterInterace::type::file);
	}

	return(this->writer->open(destination));
}

void TomahawkReader::setDetailedProgress(const bool yes){
	if(yes){
		SILENT = 0;
		this->progress.SetDetailed(true);
	} else
		this->progress.SetDetailed(false);
}

bool TomahawkReader::Open(const std::string input){
	if(input.size() == 0){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "No filename" << std::endl;
		return false;
	}

	this->stream_.open(input, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "could not open " << input << "..." << std::endl;
		return false;
	}
	this->filesize_ = this->stream_.tellg();
	this->stream_.seekg(0);

	return true;
}

bool TomahawkReader::ValidateHeader(void){
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Bad file stream..." << std::endl;
		return false;
	}

	// Read some data
	this->stream_.read(this->buffer_.data, Constants::TOMAHAWK_HEADER_LENGTH);
	this->buffer_.pointer += Constants::TOMAHAWK_HEADER_LENGTH;

	return this->Validate();
}

inline bool TomahawkReader::ValidateHeader(std::ifstream& in) const{
	char MAGIC[Constants::WRITE_HEADER_MAGIC_LENGTH];
	in.read(MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH);

	if(strncmp(MAGIC, Constants::WRITE_HEADER_MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH) == 0)
		return true;
	return false;
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
			std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Illegal seek..." << std::endl;
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
				std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Illegal seek..." << std::endl;
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
			std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Illegal seek..." << std::endl;
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
			std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Illegal seek..." << std::endl;
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
			std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Illegal seek..." << std::endl;
			return(0);
		}
		totalSize += *BSIZE;
	}

	this->stream_.seekg(this->filesize_ - Constants::eof_length*sizeof(U64) - sizeof(U32));
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Illegal seek..." << std::endl;
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

	if(!SILENT)
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

	if(!SILENT)
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

	//std::cerr << "getblock " << blockID  << " seek to " << this->totempole_[blockID].byte_offset << std::endl;
	this->stream_.seekg(this->totempole_[blockID].byte_offset);
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed search..." << std::endl;
		return false;
	}

	U32 readLength = 0;
	if(blockID != this->totempole_.getBlocks() - 1)
		readLength = this->totempole_[blockID + 1].byte_offset - this->totempole_[blockID].byte_offset;
	else
		readLength = this->filesize_ - Constants::eof_length*sizeof(U64) - this->totempole_[this->totempole_.getBlocks()-1].byte_offset;

	if(readLength > this->buffer_.capacity()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Overflowing capacity: " << readLength << '/' << this->buffer_.capacity() << std::endl;
		exit(1);
	}

	//this->manager_.Add(&this->data_.data[this->data_.pointer], this->totempole_[blockID]);
	this->blockDataOffsets_.push_back(DataOffsetPair(&this->data_.data[this->data_.pointer], this->totempole_[blockID]));
	if(!this->stream_.read(&this->buffer_.data[0], readLength)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed read: " << this->stream_.good() << '\t' << this->stream_.fail() << '/' << this->stream_.eof() << std::endl;
		//std::cerr << this->stream_.gcount() << '/' << readLength << std::endl;
		return false;
	}
	this->buffer_.pointer = readLength;

	if(!this->gzip_controller_.Inflate(this->buffer_, this->data_)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to inflate data..." << std::endl;
		return false;
	}

	//std::cerr << "Done: " << blockID << std::endl;

	return true;
}

bool TomahawkReader::Validate(void){
	if(this->buffer_.size() < Constants::WRITE_HEADER_MAGIC_LENGTH + sizeof(float) + sizeof(U64)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to validate Tomahawk header..." << std::endl;
		return false;
	}

	if(strncmp(this->buffer_.data, Constants::WRITE_HEADER_MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH) != 0){
		std::cerr << std::string(&this->buffer_.data[0], Constants::WRITE_HEADER_MAGIC_LENGTH);
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Illegal tomahawk file: failed MAGIC..." << std::endl;
		return false;
	}

	const float* version = reinterpret_cast<const float*>(&this->buffer_[Constants::WRITE_HEADER_MAGIC_LENGTH]);
	const U64* samples = reinterpret_cast<const U64*>(&this->buffer_[Constants::WRITE_HEADER_MAGIC_LENGTH + sizeof(float)]);

	this->version = *version;
	this->samples = *samples;
	this->progress.SetSamples(this->samples); // push to progressbar

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

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "BGZF") << "Inflating " << this->totempole_.getBlocks() << " blocks..." << std::endl;

	for(U32 i = 0; i < this->totempole_.getBlocks(); ++i)
		(*this.*func__)(i);

	return true;
}


} /* namespace Tomahawk */
