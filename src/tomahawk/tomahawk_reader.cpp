#include <tomahawk/tomahawk_reader.h>
#include "genotype_container.h"

namespace tomahawk {

// Remember to resize buffers to header.getLargestBlockSize()+64 after header is loaded
TomahawkReader::TomahawkReader() :
	filesize_(0),
	offset_end_of_data_(0),
	bit_width_(0),
	dropGenotypes(false),
	showHeader(true),
	index_(nullptr),
	writer(nullptr),
	interval_tree(nullptr),
	interval_tree_entries(nullptr)
{}

TomahawkReader::~TomahawkReader(){
	this->buffer_.deleteAll();
	this->data_.deleteAll();
	this->outputBuffer_.deleteAll();
	delete this->writer;
	delete this->index_;

	// Todo: fix
	//for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i)
	//	delete this->interval_tree[i];
	delete [] this->interval_tree;
	delete [] this->interval_tree_entries;
}

bool TomahawkReader::open(const std::string input){
	if(input.size() == 0){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "No input filename..." << std::endl;
		return false;
	}

	this->stream_.open(input, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream_.good()){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to open file handle: " << input << std::endl;
	}
	this->filesize_ = this->stream_.tellg();

	this->stream_.seekg(this->filesize_ - TWK_FOOTER_LENGTH);
	this->stream_ >> this->footer_;
	if(!this->stream_.good()){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Stream corrupted after loading footer..." << std::endl;
		return false;
	}

	if(this->footer_.validate() == false){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to validate footer..." << std::endl;
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
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to load header data..." << std::endl;
		return false;
	}

	if(!this->stream_.good()){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Stream is bad..." << std::endl;
		return false;
	}

	if(this->header_.validate() == false){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to validate header..." << std::endl;
		return false;
	}

	if(this->header_.magic_.major_version == 0 && this->header_.magic_.minor_version < 5){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Legacy file not supported..." << std::endl;
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
	//	std::cerr << helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << this->totempole_.getHeader().getNumberBlocks() << " blocks into " << helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

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
		std::cerr << helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << blocks.size() << " blocks into " << helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

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
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Uncompressed size is 0..." << std::endl;
		return false;
	}

	U64 totalBlocks = 0;
	for(U32 i = 0; i < blocks.size(); ++i)
		totalBlocks += blocks[i].second - blocks[i].first;

	if(!SILENT)
		std::cerr << helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << totalBlocks << " blocks into " << helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

	this->data_.resize(buffer_size + 1000);

	for(U32 i = 0; i < blocks.size(); ++i){
		for(U32 j = blocks[i].first; j < blocks[i].second; ++j){
			if(!this->getBlock(j)){
				std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Could not get block: " << j << '/' << blocks[i].second << std::endl;
				return false;
			}
		}
	}
	return true;
}


bool TomahawkReader::getBlock(const U32 blockID){
	if(!this->stream_.good()){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Bad file stream for block " << blockID << "..." << std::endl;
		return false;
	}

	this->stream_.seekg(this->index_->getContainer()[blockID].byte_offset);
	if(!this->stream_.good()){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed search..." << std::endl;
		return false;
	}

	const U32 readLength = this->index_->getContainer()[blockID].byte_offset_end - this->index_->getContainer()[blockID].byte_offset;

	if(readLength > this->buffer_.capacity()){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Overflowing capacity: " << readLength << '/' << this->buffer_.capacity() << std::endl;
		exit(1);
	}

	this->blockDataOffsets_.push_back(DataOffsetPair(&this->data_[this->data_.size()], this->index_->getContainer()[blockID].uncompressed_size, this->index_->getContainer()[blockID]));
	if(!this->stream_.read(this->buffer_.data(), readLength)){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed read: " << this->stream_.good() << '\t' << this->stream_.fail() << '/' << this->stream_.eof() << std::endl;
		//std::cerr << this->stream_.gcount() << '/' << readLength << std::endl;
		return false;
	}
	this->buffer_.n_chars = readLength;

	if(!this->tgzf_controller_.Inflate(this->buffer_, this->data_)){
		std::cerr << helpers::timestamp("ERROR", "TGZF") << "Failed to inflate data..." << std::endl;
		return false;
	}

	return true;
}

void TomahawkReader::DetermineBitWidth(void){
	if(this->header_.magic_.getNumberSamples() <= constants::UPPER_LIMIT_SAMPLES_8B - 1){
		if(!SILENT){
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " < " << constants::UPPER_LIMIT_SAMPLES_8B << "..." << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Using 8-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(BYTE);
	} else if(this->header_.magic_.getNumberSamples() <= constants::UPPER_LIMIT_SAMPLES_16B - 1){
		if(!SILENT){
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " < " << constants::UPPER_LIMIT_SAMPLES_16B << "..." << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Using 16-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(U16);
	} else if(this->header_.magic_.getNumberSamples() <= constants::UPPER_LIMIT_SAMPLES_32B - 1){
		if(!SILENT){
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " < " << constants::UPPER_LIMIT_SAMPLES_32B << "..." << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Using 32-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(U32);
	} else {
		if(!SILENT){
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " > " << constants::UPPER_LIMIT_SAMPLES_32B << "... Skip" << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->header_.magic_.getNumberSamples() << " < " << constants::UPPER_LIMIT_SAMPLES_64B << "..." << std::endl;
			std::cerr << helpers::timestamp("LOG", "RLE") << "Using 64-bit width..." << std::endl;
		}
		this->bit_width_ = sizeof(U64);
	}
}

bool TomahawkReader::outputBlocks(std::vector<U32>& blocks){
	typedef bool (tomahawk::TomahawkReader::*blockFunction)(const U32 blockID);
	blockFunction func__ = nullptr;

	switch(this->bit_width_){
	case 1: func__ = &TomahawkReader::outputBlock<BYTE>; break;
	case 2: func__ = &TomahawkReader::outputBlock<U16>;  break;
	case 4: func__ = &TomahawkReader::outputBlock<U32>;  break;
	case 8: func__ = &TomahawkReader::outputBlock<U64>;  break;
	default: exit(1); break;
	}

	if(!SILENT)
		std::cerr << helpers::timestamp("LOG","TGZF") << "Inflating " << blocks.size() << " blocks..." << std::endl;

	// Output header
	if(this->showHeader)
		std::cout << this->header_.getLiterals() + "\n##tomahawk_viewCommand=" + helpers::program_string() << std::endl;


	for(U32 i = 0; i < blocks.size(); ++i)
		(*this.*func__)(blocks[i]);

	return true;
}

bool TomahawkReader::summaryIndividuals(){
	typedef bool (tomahawk::TomahawkReader::*blockFunction)(std::vector<temp>& stats, const char* const data, const U32 blockID);
	blockFunction func__ = nullptr;

	switch(this->bit_width_){
	case 1: func__ = &TomahawkReader::statsIndividual<BYTE>; break;
	case 2: func__ = &TomahawkReader::statsIndividual<U16>;  break;
	case 4: func__ = &TomahawkReader::statsIndividual<U32>;  break;
	case 8: func__ = &TomahawkReader::statsIndividual<U64>;  break;
	default:
		std::cerr << helpers::timestamp("ERROR") << "Word sizing could not be determined!" << std::endl;
		exit(1);
		break;
	}

	std::vector<temp> stats(this->getHeader().getMagic().n_samples);
	U32 previousContigID = 0;
	for(U32 i = 0; i < this->index_->getContainer().size(); ++i){
		if(!this->getBlock(i)){
			std::cerr << "failed to get block" << std::endl;
			return false;
		}

		if(this->index_->getContainer()[i].contigID != previousContigID && i != 0){
			for(U32 i = 0; i < stats.size(); ++i){
				std::cout << previousContigID << "\t" << i << "\t" << stats[i] << std::endl;
				stats[i].clear();
			}

			previousContigID = this->index_->getContainer()[i].contigID;
		}

		(*this.*func__)(stats, this->data_.data(),i);

		// Reset buffers
		this->outputBuffer_.reset(); // reset
		this->data_.reset(); // reset

		//std::cerr << i << "/" << this->index_->getContainer().size() << std::endl;
	}

	for(U32 i = 0; i < stats.size(); ++i){
		std::cout << previousContigID << "\t" << i << "\t" << stats[i] << std::endl;
		stats[i].clear();
	}

	return true;
}

bool TomahawkReader::outputBlocks(){
	typedef bool (tomahawk::TomahawkReader::*blockFunction)(const U32 blockID);
	blockFunction func__ = nullptr;

	if(this->interval_tree_entries != nullptr){
		switch(this->bit_width_){
		case 1: func__ = &TomahawkReader::outputBlockFilter<BYTE>; break;
		case 2: func__ = &TomahawkReader::outputBlockFilter<U16>;  break;
		case 4: func__ = &TomahawkReader::outputBlockFilter<U32>;  break;
		case 8: func__ = &TomahawkReader::outputBlockFilter<U64>;  break;
		default:
			std::cerr << helpers::timestamp("ERROR") << "Word sizing could not be determined!" << std::endl;
			exit(1);
			break;
		}
	} else {
		switch(this->bit_width_){
		case 1: func__ = &TomahawkReader::outputBlock<BYTE>; break;
		case 2: func__ = &TomahawkReader::outputBlock<U16>;  break;
		case 4: func__ = &TomahawkReader::outputBlock<U32>;  break;
		case 8: func__ = &TomahawkReader::outputBlock<U64>;  break;
		default:
			std::cerr << helpers::timestamp("ERROR") << "Word sizing could not be determined!" << std::endl;
			exit(1);
			break;
		}
	}

	if(!SILENT)
		std::cerr << helpers::timestamp("LOG", "TGZF") << "Inflating " << this->index_->getContainer().size() << " blocks..." << std::endl;

	// Output header
	if(this->showHeader){
		std::cout << this->header_.getLiterals() << std::endl;
		std::cout << "##INFO=<ID=HWE_P,Number=1,Type=Float,Description=\"Hardy-Weinberg P-value (Fisher's exact test)\">" << std::endl;
		std::cout << "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">" << std::endl;
		std::cout << "##tomahawk_viewCommand=" + helpers::program_string() << std::endl;
		if(this->dropGenotypes){
			std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
		} else {
			std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";

			std::cout << this->header_.getSample(0);
			for(U32 i = 1; i < this->header_.magic_.getNumberSamples(); ++i)
				std::cout << '\t' << this->header_.getSample(i);
		}

		std::cout.put('\n');
	}

	if(this->interval_tree_entries != nullptr){
		// [a, b] overlaps with [x, y] iff b > x and a < y.
		for(U32 i = 0; i < this->index_->size(); ++i){
			// No matches in this contig
			if(this->interval_tree_entries[this->index_->getContainer().at(i).contigID].size() == 0)
				continue;

			// Find overlapping matches
			std::vector<interval_type> ret = this->interval_tree[this->index_->getContainer().at(i).contigID]->findOverlapping(this->index_->getContainer().at(i).min_position, this->index_->getContainer().at(i).max_position);

			// If there are any then use
			if(ret.size()) (*this.*func__)(i);
		}
	} else { // no filter
		for(U32 i = 0; i < this->index_->getContainer().size(); ++i){
			(*this.*func__)(i);
		}
	}

	return true;
}

bool TomahawkReader::addRegions(const std::vector<std::string>& intervals){
	// No interval strings given
	if(intervals.size() == 0)
		return true;

	// Build interval tree and interval vector is not set
	if(this->interval_tree_entries == nullptr)
		this->interval_tree_entries = new std::vector<interval_type>[this->getHeader().getMagic().getNumberContigs()];

	if(this->interval_tree == nullptr){
		this->interval_tree = new tree_type*[this->getHeader().getMagic().getNumberContigs()];
		for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i)
			this->interval_tree[i] = nullptr;
	} else delete [] this->interval_tree;

	// Parse intervals
	for(U32 i = 0; i < intervals.size(); ++i){
		interval_type interval;
		// has colon (:)
		if(intervals[i].find(':') != std::string::npos){
			std::vector<std::string> first = helpers::split(intervals[i], ':');
			if(first.size() != 2){
				std::cerr << helpers::timestamp("ERROR") << "Malformed interval string: " << intervals[i] << "..." << std::endl;
				return false;
			}

			const S32 contigID = this->header_.getContigID(first[0]);
			if(contigID == -1){
				std::cerr << helpers::timestamp("ERROR") << "Contig: " << intervals[i] << " is not defined in this file..." << std::endl;
				return false;
			}

			interval.state = algorithm::ContigInterval::INTERVAL_FULL;
			interval.contigID = contigID;
			std::vector<std::string> sections = helpers::split(first[1],'-');
			if(sections.size() == 2){ // Is contig + position->position
				if(std::regex_match(sections[0], std::regex("^[0-9]{1,}([\\.]{1}[0-9]{1,})?([eE]{1}[0-9]{1,})?$"))){
					interval.start = atof(sections[0].data());
				} else {
					std::cerr << helpers::timestamp("ERROR") << "Illegal number: " << sections[0] << " in interval string " << intervals[i] << std::endl;
					return false;
				}

				if(std::regex_match(sections[1], std::regex("^[0-9]{1,}([\\.]{1}[0-9]{1,})?([eE]{1}[0-9]{1,})?$"))){
					interval.stop = atof(sections[1].data());
				} else {
					std::cerr << helpers::timestamp("ERROR") << "Illegal number: " << sections[1] << " in interval string " << intervals[i] << std::endl;
					return false;
				}

			} else if(sections.size() == 1){ // Is contig + position
				interval.state = algorithm::ContigInterval::INTERVAL_POSITION;
				if(std::regex_match(sections[0], std::regex("^[0-9]{1,}([\\.]{1}[0-9]{1,})?([eE]{1}[0-9]{1,})?$"))){
					interval.start = atof(sections[0].data());
					interval.stop  = atof(sections[0].data());
				} else {
					std::cerr << helpers::timestamp("ERROR") << "Illegal number: " << sections[0] << " in interval string " << intervals[i] << std::endl;
					return false;
				}

			} else {
				std::cerr << helpers::timestamp("ERROR") << "Malformed interval string: " << intervals[i] << std::endl;
				return false;
			}
		} else { // Is contig  only
			interval.state = algorithm::ContigInterval::INTERVAL_CONTIG_ONLY;
			const S32 contigID = this->header_.getContigID(intervals[i]);
			if(contigID == -1){
				std::cerr << helpers::timestamp("ERROR") << "Contig: " << intervals[i] << " is not defined in this file..." << std::endl;
				return false;
			}
			interval.contigID = contigID;
			interval.start = 0;
			interval.stop = this->header_.contigs_[contigID].n_bases + 1;
		}
		// Store interval
		this->interval_tree_entries[interval.contigID].push_back(interval_type(interval));
	}

	for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i){
		delete this->interval_tree[i];

		if(this->interval_tree_entries[i].size() != 0){
			this->interval_tree[i] = new tree_type(this->interval_tree_entries[i]);
		} else
			this->interval_tree[i] = nullptr;
	}

	return true;
}


} /* namespace Tomahawk */
