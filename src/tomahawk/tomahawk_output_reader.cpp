#include <cstdlib>
#include <cstring>
#include <iostream>

#include "io/output_writer.h"
#include "tomahawk/tomahawk_output_reader.h"
#include "io/compression/gz_constants.h"
#include "io/compression/gz_header.h"
#include "support/helpers.h"
#include "math/output_statistics.h"
#include "tomahawk_output_stats.h"
#include "aggregation_parameters.h"

namespace tomahawk {

TomahawkOutputReader::TomahawkOutputReader() :
	filesize_(0),
	offset_end_of_data_(0),
	index_(nullptr),
	buffer_(3000000),
	data_(3000000),
	outputBuffer_(3000000),
	interval_tree(nullptr),
	interval_tree_entries(nullptr),
	writer_(nullptr)
{

}

TomahawkOutputReader::~TomahawkOutputReader(){
	delete this->index_;
	delete [] this->interval_tree_entries;
	if(this->interval_tree != nullptr){
		for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i){
			delete this->interval_tree[i];
		}
		delete [] this->interval_tree;
	}
	delete this->writer_;
}

bool TomahawkOutputReader::open(const std::string input){
	if(input.size() == 0){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "No input filename..." << std::endl;
		return false;
	}
	this->parameters_.input_file = input;

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
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to validate footer! The file is truncated or corrupted..." << std::endl;
		return false;
	}

	// Seek to start of index
	this->stream_.seekg(this->footer_.offset_end_of_data);
	const U32 l_index_data = (this->filesize_ - TWK_FOOTER_LENGTH) - this->stream_.tellg();
	buffer_type index_buffer(l_index_data + 1024);
	this->stream_.read(index_buffer.data(), l_index_data);
	index_buffer.n_chars = l_index_data;
	this->index_ = new index_type(index_buffer.data(), index_buffer.size());

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

	if(this->header_.magic_.major_version == 0 && this->header_.magic_.minor_version < 0.4){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Legacy file not supported..." << std::endl;
		return false;
	}

	this->offset_end_of_data_ = this->footer_.offset_end_of_data;

	return true;
}

bool TomahawkOutputReader::printHeader(std::ostream& stream, std::vector<std::string>& extra){
	for(U32 i = 0; i < extra.size(); ++i){
		this->getHeader().getLiterals() += "\n" + extra[i];
	}

	return(this->printHeader(stream));
}

bool TomahawkOutputReader::printHeader(std::ostream& stream) const{
	if(this->parameters_.showHeader == true){
		std::cout << this->getHeader().getLiterals() << '\n';
		std::cout << "FLAG\tCHROM_A\tPOS_A\tCHROM_B\tPOS_B\tREF_REF\tREF_ALT\tALT_REF\tALT_ALT\tD\tDprime\tR\tR2\tP\tChiSqModel\tChiSqTable\n";
	}

	return true;
}

int TomahawkOutputReader::nextBlock(const bool clear, const bool clear_raw){
	// Stream died
	if(this->stream_.good() == false){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Stream died!" << std::endl;
		return -1;
	}

	// EOF
	// tellg will always return a positive value here
	// or it would've failed at good() check
	if((U64)this->stream_.tellg() == this->offset_end_of_data_)
		return 0;

	// Reset compressed_buffer
	if(clear_raw) this->buffer_.reset();

	// Read TGZF header
	this->buffer_.resize(sizeof(tgzf_header_type));
	this->stream_.read(this->buffer_.data(),  io::constants::TGZF_BLOCK_HEADER_LENGTH);
	const tgzf_header_type* h = reinterpret_cast<const tgzf_header_type*>(this->buffer_.data());
	this->buffer_.n_chars = io::constants::TGZF_BLOCK_HEADER_LENGTH;
	if(!h->Validate()){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Failed to validate!" << std::endl;
		return -2;
	}

	this->buffer_.resize(h->BSIZE); // make sure all data will fit

	// Recast because if compressed_buffer is (actually) resized then the pointer address is incorrect
	// resulting in segfault
	h = reinterpret_cast<const tgzf_header_type*>(this->buffer_.data());

	this->stream_.read(&this->buffer_.buffer[io::constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - io::constants::TGZF_BLOCK_HEADER_LENGTH);
	if(!this->stream_.good()){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Truncated file..." << std::endl;
		return -3;
	}

	this->buffer_.n_chars = h->BSIZE;
	const U32 uncompressed_size = *reinterpret_cast<const U32*>(&this->buffer_[this->buffer_.size() - sizeof(U32)]);

	// Clear output compressed_buffer
	if(clear) {
		this->data_.reset();
		this->data_.resize(uncompressed_size);
	} else { // Otherwise resize to permit data
		this->data_.resize(this->data_.size() + uncompressed_size);
	}

	if(!this->tgzf_controller_.Inflate(this->buffer_, this->data_)){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Failed inflate!" << std::endl;
		return -4;
	}

	if(this->data_.size() == 0){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Empty data!" << std::endl;
		return 0;
	}

	// Reset iterator position and size
	//this->iterator_position_block = 0;
	//this->size = this->data_.size() / sizeof(entry_type);

	// Validity check
	if(this->data_.size() % sizeof(entry_type) != 0){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Data is corrupted!" << std::endl;
		return -5;
	}

	return 1;
}

int TomahawkOutputReader::parseBlock(std::ifstream& stream,
                                       buffer_type& inflate_buffer,
                                       buffer_type& data_buffer,
                              tgzf_controller_type& compression_manager,
                                        const bool  clear) const
{
	// Stream died
	if(stream.good() == false){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Stream died!" << std::endl;
		return -1;
	}

	// EOF
	// tellg will always return a positive value here
	// or it would've failed at good() check
	if((U64)stream.tellg() == this->offset_end_of_data_)
		return 0;

	// Read TGZF header
	inflate_buffer.resize(sizeof(tgzf_header_type));
	stream.read(inflate_buffer.data(),  io::constants::TGZF_BLOCK_HEADER_LENGTH);
	const tgzf_header_type* h = reinterpret_cast<const tgzf_header_type*>(inflate_buffer.data());
	inflate_buffer.n_chars = io::constants::TGZF_BLOCK_HEADER_LENGTH;
	if(!h->Validate()){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Failed to validate!" << std::endl;
		return -2;
	}

	inflate_buffer.resize(h->BSIZE); // make sure all data will fit

	// Recast because if compressed_buffer is (actually) resized then the pointer address is incorrect
	// resulting in segfault
	h = reinterpret_cast<const tgzf_header_type*>(inflate_buffer.data());

	stream.read(&inflate_buffer.buffer[io::constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - io::constants::TGZF_BLOCK_HEADER_LENGTH);
	if(!stream.good()){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Truncated file..." << std::endl;
		return -3;
	}

	inflate_buffer.n_chars = h->BSIZE;
	const U32 uncompressed_size = *reinterpret_cast<const U32*>(&inflate_buffer[inflate_buffer.size() - sizeof(U32)]);

	// Clear output compressed_buffer
	if(clear) {
		data_buffer.reset();
		data_buffer.resize(uncompressed_size);
	} else { // Otherwise resize to permit data
		data_buffer.resize(data_buffer.size() + uncompressed_size);
	}

	if(!this->tgzf_controller_.Inflate(inflate_buffer, data_buffer)){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Failed inflate!" << std::endl;
		return -4;
	}

	if(data_buffer.size() == 0){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Empty data!" << std::endl;
		return 0;
	}

	// Reset compressed_buffer
	inflate_buffer.reset();

	// Reset iterator position and size
	//this->iterator_position_block = 0;
	//this->size = data_buffer.size() / sizeof(entry_type);

	// Validity check
	if(data_buffer.size() % sizeof(entry_type) != 0){
		std::cerr << helpers::timestamp("ERROR", "TWO") << "Data is corrupted!" << std::endl;
		return -5;
	}

	return 1;
}

containers::OutputContainer TomahawkOutputReader::getContainerVariants(const U64 n_variants){
	size_t n_variants_loaded = 0;
	this->data_.reset();
	this->data_.resize(n_variants*sizeof(entry_type) + 65536); // make room for data
	while(true){
		if(!this->nextBlock(false))
			break;

		n_variants_loaded = this->data_.size() / sizeof(entry_type);
		//std::cerr << n_variants_loaded << "/" << n_variants << '\t' << this->data_.size() << std::endl;
		if(n_variants_loaded >= n_variants)
			break;
	}

	return(containers::OutputContainer(this->data_));
}

containers::OutputContainer TomahawkOutputReader::getContainerBytes(const size_t l_data){
	const U64 start_position = this->stream_.tellg();
	this->data_.reset();
	this->data_.resize(l_data + 65536); // make room for data
	U64 data_loaded = 0;
	while(true){
		if(!this->nextBlock(false))
			break;

		data_loaded = (U64)this->stream_.tellg() - start_position;
		if(data_loaded >= l_data)
			break;

	}

	return(containers::OutputContainer(this->data_));
}

containers::OutputContainer TomahawkOutputReader::getContainerBlocks(const U32 n_blocks){
	this->data_.reset();
	for(U32 i = 0; i < n_blocks; ++i){
		if(!this->nextBlock(false))
			break;
	}

	return(containers::OutputContainer(this->data_));
}

bool TomahawkOutputReader::seekBlock(const U32 blockID){
	if(blockID > this->getIndex().getContainer().size()){
		std::cerr << helpers::timestamp("ERROR","TOI") << "Illegal blockID (" << blockID << ">" << this->getIndex().getContainer().size() << ")!" << std::endl;
		return false;
	}

	if(!this->stream_.good()){
		std::cerr << helpers::timestamp("ERROR","TWO") << "Stream is bad!" << std::endl;
		return false;
	}

	this->stream_.seekg(this->getIndex().getContainer()[blockID].byte_offset);
	if(!this->stream_.good()){
		std::cerr << helpers::timestamp("ERROR","TWO") << "Stream is bad following seek!" << std::endl;
		return false;
	}

	return(true);
}

bool TomahawkOutputReader::seekBlock(std::ifstream& stream, const U32 blockID) const{
	if(blockID > this->getIndex().getContainer().size()){
		std::cerr << helpers::timestamp("ERROR","TOI") << "Illegal blockID (" << blockID << ">" << this->getIndex().getContainer().size() << ")!" << std::endl;
		return false;
	}

	if(!stream.good()){
		std::cerr << helpers::timestamp("ERROR","TWO") << "Stream is bad!" << std::endl;
		return false;
	}

	stream.seekg(this->getIndex().getContainer()[blockID].byte_offset);
	if(!stream.good()){
		std::cerr << helpers::timestamp("ERROR","TWO") << "Stream is bad following seek!" << std::endl;
		return false;
	}

	return(true);
}


containers::OutputContainerReference TomahawkOutputReader::getContainerReferenceBlock(const U32 blockID){
	if(!this->seekBlock(blockID)){
		this->nextBlock();
		return(containers::OutputContainerReference());
	} else {
		this->nextBlock();
		return(containers::OutputContainerReference(this->data_));
	}
}

containers::OutputContainer TomahawkOutputReader::getContainerBlock(const U32 blockID){
	if(!this->seekBlock(blockID)){
		this->nextBlock();
		return(containers::OutputContainer());
	} else {
		this->nextBlock();
		return(containers::OutputContainer(this->data_));
	}
}

containers::OutputContainerReference TomahawkOutputReader::getContainerReferenceBlock(std::vector<U32> blocks){
	if(!this->seekBlock(blocks[0])){
		this->nextBlock();
		return(containers::OutputContainerReference());
	} else {
		for(U32 i = 0; i < blocks.size(); ++i){
			if(!this->nextBlock(false))
				break;
		}
		return(containers::OutputContainerReference(this->data_));
	}
}

containers::OutputContainer TomahawkOutputReader::getContainerBlock(std::vector<U32> blocks){
	if(!this->seekBlock(blocks[0])){
		this->nextBlock();
		return(containers::OutputContainer());
	} else {
		for(U32 i = 0; i < blocks.size(); ++i){
			if(!this->nextBlock(false))
				break;
		}
		return(containers::OutputContainer(this->data_));
	}
}

bool TomahawkOutputReader::addRegions(std::vector<std::string>& positions){
	if(positions.size() == 0)
		return true;

	// No contigs (usually means data has not been loaded)
	if(this->getHeader().getMagic().getNumberContigs() == 0)
		return false;

	// Construct interval tree and interval vector if not set
	if(this->interval_tree_entries == nullptr){
		this->interval_tree_entries = new std::vector<interval_type>[this->getHeader().getMagic().getNumberContigs()];

		// Reserve some memory
		// Linked reads require that he pointers between entries do not change
		for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i)
			this->interval_tree_entries[i].reserve(1000);
	}

	if(this->interval_tree == nullptr){
		this->interval_tree = new tree_type*[this->getHeader().getMagic().getNumberContigs()];
		for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i)
			this->interval_tree[i] = nullptr;
	}

	// Parse and add region
	if(!this->__addRegions(positions))
		return false;

	for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i){
		if(this->interval_tree_entries[i].size() != 0){
			this->interval_tree[i] = new tree_type(this->interval_tree_entries[i]);
		} else
			this->interval_tree[i] = nullptr;
	}

	return true;
}

bool TomahawkOutputReader::__addRegions(std::vector<std::string>& positions){
	// For each potential interval string in vector
	for(U32 i = 0; i < positions.size(); ++i){
		if(positions[i].find(',') != std::string::npos){ // Test if string has a comma
			std::vector<std::string> ret = helpers::split(positions[i], ',');
			if(ret.size() == 1){ //
				std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
				return false;

			} else if(ret.size() == 2){
				// parse left
				interval_type intervalLeft;
				interval_type intervalRight;

				if(this->__ParseRegion(ret[0], intervalLeft) == false){
					std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Failed interpret left interval in pair!" << std::endl;
					return false;
				}

				// parse right
				if(this->__ParseRegion(ret[1], intervalRight) == false){
					std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Failed interpret right interval in pair!" << std::endl;
					return false;
				}

				// Assert that left is smaller than right
				if(intervalRight < intervalLeft)
					std::swap(intervalRight, intervalLeft);


				// Link intervals together
				this->interval_tree_entries[intervalLeft.contigID].push_back(interval_type(intervalLeft));
				interval_type& left_pointer = this->interval_tree_entries[intervalLeft.contigID].back();
				this->interval_tree_entries[intervalRight.contigID].push_back(interval_type(intervalRight));
				interval_type& right_pointer = this->interval_tree_entries[intervalRight.contigID].back();
				right_pointer.value = &left_pointer;
				left_pointer.value  = &right_pointer;

			} else {
				std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
				return false;
			}
		}
		// Has no comma in string
		else {
			interval_type interval;
			if(this->__ParseRegion(positions[i], interval))
				this->interval_tree_entries[interval.contigID].push_back(interval_type(interval));
		}
	}

	return true;
}

bool TomahawkOutputReader::__ParseRegion(const std::string& region, interval_type& interval) const{
	if(region.size() == 0)
		return false;

	// Search for colon
	std::vector<std::string> ret = helpers::split(region, ':');

	// If vector does not contain a colon
	if(ret.size() == 1){
		if(ret[0].find('-') != std::string::npos){
			std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}

		// is contigID only
		const S32 contigID = this->getHeader().getContigID(ret[0]);
		if(contigID < 0){
			std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << region << " is not defined in the header!" << std::endl;
			return false;
		}

		interval(contigID, 0, this->getHeader().contigs_[contigID].n_bases);
		interval.state = interval_type::INTERVAL_TYPE::INTERVAL_CONTIG_ONLY;
	}
	// If vector contain colon
	else if(ret.size() == 2){
		// is contigID:pos-pos
		const S32 contigID = this->getHeader().getContigID(ret[0]);
		if(contigID < 0){
			std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << ret[0] << " is not defined in the header!" << std::endl;
			return false;
		}

		std::vector<std::string> retPos = helpers::split(ret[1], '-');
		if(retPos.size() == 1){
			// only one pos
			const double pos = std::stod(retPos[0]);
			//std::cerr << "single position: " << pos << std::endl;
			interval(contigID, pos, pos);
			interval.state = interval_type::INTERVAL_TYPE::INTERVAL_POSITION;

		} else if(retPos.size() == 2){
			// is two positions
			double posA = std::stod(retPos[0]);
			double posB = std::stod(retPos[1]);

			// Swap pA and pB iff pB > pA
			if(posB < posA)
				std::swap(posA, posB);

			interval(contigID, posA, posB);
			interval.state = interval_type::INTERVAL_TYPE::INTERVAL_FULL;
		} else {
			std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}
	}
	// contains > 1 colons
	// illegal
	else {
		std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::view(void){
	if(!this->openWriter()) return false;

	if(this->interval_tree != nullptr) // If regions have been set: use region-filter function
		return(this->__viewRegion());
	else if(this->filters_.any_filter_user_set)
		return(this->__viewFilter()); // Otherwise normal filter function
	 else
		return(this->__viewOnly());
}

bool TomahawkOutputReader::__viewOnly(void){
	std::vector<std::string> extra;
	extra.push_back("##tomahawk_viewCommand=" + helpers::program_string());
	extra.push_back("##tomahawk_viewFilters=" + this->filters_.getInterpretedString() + " filter=NO regions=NO");

	//if(this->showHeader_ == true)
	//	this->printHeader(std::cout, extra);

	const std::string version_string = std::to_string(this->header_.magic_.major_version) + "." + std::to_string(this->header_.magic_.minor_version) + "." + std::to_string(this->header_.magic_.patch_version);

	assert(this->writer_ != nullptr);
	this->writer_->writeHeaders(this->getHeader());

	// Natural output required parsing
	typedef std::ostream& (entry_type::*func)(std::ostream& os, const contig_type* const contigs) const;

	func ld_write_function = &entry_type::write;

	//this->parseBlock();
	if(this->parameters_.output_type == TWK_OUTPUT_LD){
		while(this->nextBlock()){
			containers::OutputContainerReference o = this->getContainerReference();
			(o[0].*ld_write_function)(std::cout, this->getHeader().contigs_);

			for(U32 i = 1; i < o.size(); ++i)
				(o[i].*ld_write_function)(std::cout, this->getHeader().contigs_);

		}
	} else {
		while(this->nextBlock()){
			containers::OutputContainerReference o = this->getContainerReference();
			*this->writer_ << o[0];

			for(U32 i = 1; i < o.size(); ++i)
				*this->writer_ << o[i];
		}
	}

	this->writer_->flush();
	if(this->writer_->isSorted())
		this->writer_->getIndex()->buildMetaIndex(this->getHeader().magic_.n_contigs);

	this->writer_->writeFinal();
	this->writer_->flush();

	return true;
}

bool TomahawkOutputReader::__viewRegion(void){
	std::vector<std::string> extra;
	extra.push_back("##tomahawk_viewCommand=" + helpers::program_string());
	if(this->filters_.any_filter_user_set)
		extra.push_back("##tomahawk_viewFilters=" + this->filters_.getInterpretedString() + " filter=YES regions=YES");

	//if(this->parameters_.showHeader == true)
	//	this->printHeader(std::cout, extra);

	if(this->interval_tree == nullptr)
		return false;

	assert(this->writer_ != nullptr);
	this->writer_->writeHeaders(this->getHeader());

	typedef bool (TomahawkOutputReader::*func_slice)(const entry_type& entry);
	func_slice region_function = &TomahawkOutputReader::__checkRegionUnsorted;
	if(this->getIndex().isSorted()) region_function = &TomahawkOutputReader::__checkRegionSorted;

	if(this->getIndex().getController().isSorted){
		if(!SILENT)
			std::cerr << helpers::timestamp("LOG") << "Using sorted query..." << std::endl;

		region_function = &TomahawkOutputReader::__checkRegionSorted;

		for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i){ // foreach contig
			//std::cerr << "i: " << this->interval_tree_entries[i].size() << std::endl;
			for(U32 j = 0; j < this->interval_tree_entries[i].size(); ++j){
				//std::cerr << this->interval_tree_entries[i][j] << std::endl;
				std::vector<U32> blocks = this->index_->findOverlaps(this->interval_tree_entries[i][j].contigID, this->interval_tree_entries[i][j].start, this->interval_tree_entries[i][j].stop);
				for(U32 b = 0; b < blocks.size(); ++b){
					output_container_reference_type o = this->getContainerReferenceBlock(blocks[b]);
					for(U32 i = 0; i < o.size(); ++i)
						(this->*region_function)(o[i]);

				} // end of blocks
			} // end of intervals in contig
		} // end of contigs
	} // end case sorted
	else {
		if(!SILENT)
			std::cerr << helpers::timestamp("LOG") << "Using unsorted query..." << std::endl;

		while(this->nextBlock()){
			output_container_reference_type o(this->data_);
			for(U32 i = 0; i < o.size(); ++i)
				(this->*region_function)(o[i]);
		}
	}

	this->writer_->flush();
	if(this->writer_->isSorted())
		this->writer_->getIndex()->buildMetaIndex(this->getHeader().magic_.n_contigs);
	this->writer_->writeFinal();
	this->writer_->flush();

	return true;
}

bool TomahawkOutputReader::__viewFilter(void){
	std::vector<std::string> extra;
	extra.push_back("##tomahawk_viewCommand=" + helpers::program_string());
	extra.push_back("##tomahawk_viewFilters=" + this->filters_.getInterpretedString() + " filter=YES regions=NO");

	//if(this->parameters_.showHeader == true)
	//	this->printHeader(std::cout, extra);

	assert(this->writer_ != nullptr);
	this->writer_->writeHeaders(this->getHeader());

	while(this->nextBlock()){
		output_container_reference_type o(this->data_);
		for(U32 i = 0; i < o.size(); ++i){
			if(this->filters_.filter(o[i])){
				o[i].write(std::cout, this->getHeader().contigs_);
			}
		}
	} // end while next block

	this->writer_->flush();
	if(this->writer_->isSorted())
		this->writer_->getIndex()->buildMetaIndex(this->getHeader().magic_.n_contigs);
	this->writer_->writeFinal();
	this->writer_->flush();

	return true;
}

bool TomahawkOutputReader::__checkRegionSorted(const entry_type& entry){
	typedef std::ostream& (entry_type::*func)(std::ostream& os, const contig_type* const contigs) const;
	func ld_output_function = &entry_type::write;

	// If iTree for contigA exists
	if(this->interval_tree[entry.AcontigID] != nullptr){
		std::vector<interval_type> overlaps_in_tree = this->interval_tree[entry.AcontigID]->findOverlapping(entry.Aposition, entry.Aposition);
		if(overlaps_in_tree.size() > 0){
			for(U32 i = 0; i < overlaps_in_tree.size(); ++i){
				if(overlaps_in_tree[i].value != nullptr){ // if linked
					if((entry.BcontigID == overlaps_in_tree[i].value->contigID) &&
					   (entry.Bposition >= overlaps_in_tree[i].value->start &&
						entry.Bposition <= overlaps_in_tree[i].value->stop))
					{
						if(this->filters_.filter(entry)){
							//entry.write(std::cout, this->contigs);
							//*this->writer << entry;
							if(this->parameters_.output_type == TWK_OUTPUT_LD) (entry.*ld_output_function)(std::cout, this->getHeader().contigs_);
							else *this->writer_ << entry;
						}

						return true;
					} // end match
				} else { //  not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						if(this->parameters_.output_type == TWK_OUTPUT_LD) (entry.*ld_output_function)(std::cout, this->getHeader().contigs_);
						else *this->writer_ << entry;
					}

					return true;
				}
			}
		}
	}

	return false;
}


bool TomahawkOutputReader::__checkRegionUnsorted(const entry_type& entry){
	typedef std::ostream& (entry_type::*func)(std::ostream& os, const contig_type* const contigs) const;
	func ld_output_function = &entry_type::write;

	// If iTree for contigA exists
	if(this->interval_tree[entry.AcontigID] != nullptr){
		std::vector<interval_type> overlaps_in_tree = this->interval_tree[entry.AcontigID]->findOverlapping(entry.Aposition, entry.Aposition);
		if(overlaps_in_tree.size() > 0){
			for(U32 i = 0; i < overlaps_in_tree.size(); ++i){
				if(overlaps_in_tree[i].value != nullptr){ // if linked
					if((entry.BcontigID == overlaps_in_tree[i].value->contigID) &&
					   (entry.Bposition >= overlaps_in_tree[i].value->start &&
					    entry.Bposition <= overlaps_in_tree[i].value->stop))
					{
						if(this->filters_.filter(entry)){
							//entry.write(std::cout, this->contigs);
							//*this->writer << entry;
							if(this->parameters_.output_type == TWK_OUTPUT_LD) (entry.*ld_output_function)(std::cout, this->getHeader().contigs_);
							else *this->writer_ << entry;
						}

						return true;
					} // end match
				} else { //  not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						if(this->parameters_.output_type == TWK_OUTPUT_LD) (entry.*ld_output_function)(std::cout, this->getHeader().contigs_);
						else *this->writer_ << entry;
					}

					return true;
				}
			}
		}
	}

	// If iTree for contigB exists
	if(this->interval_tree[entry.BcontigID] != nullptr){
		std::vector<interval_type> overlaps_in_tree = this->interval_tree[entry.BcontigID]->findOverlapping(entry.Bposition, entry.Bposition);
		if(overlaps_in_tree.size() > 0){
			for(U32 i = 0; i < overlaps_in_tree.size(); ++i){
				if(overlaps_in_tree[i].value != nullptr){ // if linked
					if((entry.AcontigID == overlaps_in_tree[i].value->contigID) &&
					   (entry.Aposition >= overlaps_in_tree[i].value->start &&
						entry.Aposition <= overlaps_in_tree[i].value->stop)){
						if(this->filters_.filter(entry)){
							//entry.write(std::cout, this->contigs);
							//*this->writer << entry;
							if(this->parameters_.output_type == TWK_OUTPUT_LD) (entry.*ld_output_function)(std::cout, this->getHeader().contigs_);
							else *this->writer_ << entry;
						}
						return true;
					} // end match
				} else { // not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						if(this->parameters_.output_type == TWK_OUTPUT_LD) (entry.*ld_output_function)(std::cout, this->getHeader().contigs_);
						else *this->writer_ << entry;
					}

					return true;
				}
			}
		} // end if any hit in iTree b
	} // end iTree b

	return false;
}

bool TomahawkOutputReader::__concat(const std::vector<std::string>& files, const std::string& output){
	if(files.size() == 0){
		std::cerr << helpers::timestamp("ERROR","TWO") << "No input files..." << std::endl;
		return false;
	}

	// open first one
	if(!SILENT)
		std::cerr << helpers::timestamp("LOG", "CONCAT") << "Opening input: " << files[0] << " (" << 1 << "/" << files.size() << ")" << "..." << std::endl;

	if(!this->open(files[0])){
		std::cerr << helpers::timestamp("ERROR","TWO") << "Failed to parse: " << files[0] << " (" << 1 << "/" << files.size() << ")" << "..." << std::endl;
		return false;
	}

	this->getHeader().getLiterals() += "\n##tomahawk_concatCommand=" + helpers::program_string();
	this->getHeader().getLiterals() += "\n##tomahawk_concatFiles=";
	for(U32 i = 0; i < files.size(); ++i)
		this->getHeader().getLiterals() += files[i] + ',';

	io::OutputWriterBinaryFile writer;
	if(!writer.open(output)){
		std::cerr << helpers::timestamp("ERROR","SORT") << "Failed to open: " << output << "..." << std::endl;
		return false;
	}
	writer.writeHeaders(this->getHeader());

	while(this->nextBlock()){
		writer.writePrecompressedBlock(this->buffer_, this->data_.size());
		//writer << this->data_;
	}

	for(U32 i = 1; i < files.size(); ++i){
		if(!SILENT)
			std::cerr << helpers::timestamp("LOG", "CONCAT") << "Opening input: " << files[i] << " (" << i+1 << "/" << files.size() << ")" << "..." << std::endl;

		self_type second_reader;
		if(!second_reader.open(files[i])){
			std::cerr << helpers::timestamp("ERROR","TWO") << "Failed to parse: " << files[i] << " (" << i+1 << "/" << files.size() << ")" << "..." << std::endl;
			return false;
		}

		if(!(second_reader.getHeader() == this->getHeader())){
			std::cerr << helpers::timestamp("ERROR","CONCAT") << "Header mismatch: these files appear to be originating from different TWK files..." << std::endl;
			return false;
		}

		while(second_reader.nextBlock()){
			//writer << second_reader.data_;
			writer.writePrecompressedBlock(second_reader.buffer_, second_reader.data_.size());
		}

		writer.flush();
	}

	writer.setSorted(false);
	writer.setPartialSorted(false);
	writer.flush();
	writer.writeFinal();

	return true;
}


bool TomahawkOutputReader::concat(const std::vector<std::string>& files, const std::string& output){
	if(files.size() == 0){
		std::cerr << helpers::timestamp("ERROR","TWO") << "No input files given..." << std::endl;
		return false;
	}

	return(this->__concat(files, output));
}

bool TomahawkOutputReader::concat(const std::string& file_list, const std::string& output){
	if(file_list.size() == 0){
		std::cerr << helpers::timestamp("ERROR","TWO") << "No input file list given..." << std::endl;
		return false;
	}

	std::ifstream file_list_read(file_list);
	if(!file_list_read.good()){
		std::cerr << helpers::timestamp("ERROR","TWO") << "Failed to get file list..." << std::endl;
		return false;
	}

	std::vector<std::string> files;
	std::string line;
	while(getline(file_list_read, line)){
		if(line.size() == 0){
			std::cerr << helpers::timestamp("WARNING","TWO") << "Empty line. Attempting to continue..." << std::endl;
			continue;
		}
		files.push_back(line);
	}

	return(this->__concat(files, output));
}

bool TomahawkOutputReader::statistics(void){
	if(this->getIndex().getController().isSorted == false){
		std::cerr << helpers::timestamp("LOG","STATS") << "File has to be sorted..." << std::endl;
		return false;
	}

	if(this->nextBlock() == false){
		std::cerr << helpers::timestamp("ERROR") << "No valid data!" << std::endl;
		return false;
	}

	SummaryStatisticsObject statsContig; // contig vs everything
	SummaryStatisticsObject statsContigSelf; // contig vs self only
	SummaryStatisticsObject statsContigPosition; // contig-pos vs everything
	SummaryStatisticsObject statsContigPositionContig; // contig-pos-contig vs everything

	std::vector<SummaryStatisticsObject> contig_data;
	std::vector<SummaryStatisticsObject> contig_self_data;
	std::vector<SummaryStatisticsObject> contig_position_data;
	std::vector<SummaryStatisticsObject> contig_position_contig_data;

	std::vector< std::vector< SummaryStatisticsObject > > pairwise_contigs(this->getHeader().magic_.n_contigs, std::vector< SummaryStatisticsObject >(this->getHeader().magic_.n_contigs));

	const std::string CONTIG_ALL  = "contig-all";
	const std::string CONTIG_POS  = "contig-pos";
	const std::string CONTIG_POS_CONTIG = "contig-pos-contig";
	const std::string CONTIG_SELF = "contig-self";

	// Now have first data
	containers::OutputContainerReference o = this->getContainerReference();
	U32 ref_contigA   = o[0].AcontigID;
	U32 ref_positionA = o[0].Aposition;
	U32 ref_contigB   = o[0].BcontigID;

	statsContig += o[0];
	statsContigSelf += o[0];
	statsContigPosition += o[0];
	statsContigPositionContig += o[0];
	pairwise_contigs[o[0].AcontigID][o[0].BcontigID] += o[0];

	for(U32 i = 1; i < o.size(); ++i){
		if(o[i].AcontigID != ref_contigA){
			//statsContig.print(std::cout, CONTIG_ALL, this->getHeader().contigs_[ref_contigA].name, -1, ".");
			contig_data.push_back(statsContig);
			statsContig.reset();
			//statsContigSelf.print(std::cout, CONTIG_SELF, this->getHeader().contigs_[ref_contigA].name, -1, this->getHeader().contigs_[ref_contigA].name);
			contig_self_data.push_back(statsContigSelf);
			statsContigSelf.reset();
			//statsContigPosition.print(std::cout, CONTIG_POS, this->getHeader().contigs_[ref_contigA].name, ref_positionA, ".");
			contig_position_data.push_back(statsContigPosition);
			statsContigPosition.reset();
			ref_positionA = o[i].Aposition;
			ref_contigA   = o[i].AcontigID;
			ref_contigB   = o[i].BcontigID;
		}

		if(o[i].AcontigID != ref_contigA || o[i].Aposition != ref_positionA){
			//statsContigPosition.print(std::cout, CONTIG_POS, this->getHeader().contigs_[ref_contigA].name, ref_positionA, ".");
			contig_position_data.push_back(statsContigPosition);
			statsContigPosition.reset();
			ref_positionA = o[i].Aposition;
		}

		if(o[i].AcontigID != ref_contigA || o[i].Aposition != ref_positionA || o[i].BcontigID != ref_contigB){
			//statsContigPositionContig.print(std::cout, CONTIG_POS_CONTIG, this->getHeader().contigs_[ref_contigA].name, ref_positionA, this->getHeader().contigs_[ref_contigB].name);
			contig_position_contig_data.push_back(statsContigPositionContig);
			statsContigPositionContig.reset();
			ref_contigB = o[i].BcontigID;
		}

		statsContig += o[i];
		statsContigSelf += o[i];
		statsContigPosition += o[i];
		statsContigPositionContig += o[i];
		pairwise_contigs[o[i].AcontigID][o[i].BcontigID] += o[i];
	}

	// For the remainder of the blocks
	while(this->nextBlock()){
		containers::OutputContainerReference o = this->getContainerReference();

		for(U32 i = 0; i < o.size(); ++i){
			if(o[i].AcontigID != ref_contigA){
				//statsContig.print(std::cout, CONTIG_ALL, this->getHeader().contigs_[ref_contigA].name, -1, ".");
				contig_data.push_back(statsContig);
				statsContig.reset();
				//statsContigSelf.print(std::cout, CONTIG_SELF, this->getHeader().contigs_[ref_contigA].name, -1, this->getHeader().contigs_[ref_contigA].name);
				contig_self_data.push_back(statsContigSelf);
				statsContigSelf.reset();
				//statsContigPosition.print(std::cout, CONTIG_POS, this->getHeader().contigs_[ref_contigA].name, ref_positionA, ".");
				contig_position_data.push_back(statsContigPosition);
				statsContigPosition.reset();
				ref_positionA = o[i].Aposition;
				ref_contigA   = o[i].AcontigID;
				ref_contigB   = o[i].BcontigID;
			}

			if(o[i].AcontigID != ref_contigA || o[i].Aposition != ref_positionA){
				//statsContigPosition.print(std::cout, CONTIG_POS, this->getHeader().contigs_[ref_contigA].name, ref_positionA, ".");
				contig_position_data.push_back(statsContigPosition);
				statsContigPosition.reset();
				ref_positionA = o[i].Aposition;
			}

			if(o[i].AcontigID != ref_contigA || o[i].Aposition != ref_positionA || o[i].BcontigID != ref_contigB){
				//statsContigPositionContig.print(std::cout, CONTIG_POS_CONTIG, this->getHeader().contigs_[ref_contigA].name, ref_positionA, this->getHeader().contigs_[ref_contigB].name);
				contig_position_contig_data.push_back(statsContigPositionContig);
				statsContigPositionContig.reset();
				ref_contigB = o[i].BcontigID;
			}

			statsContig += o[i];
			statsContigSelf += o[i];
			statsContigPosition += o[i];
			statsContigPositionContig += o[i];
			pairwise_contigs[o[i].AcontigID][o[i].BcontigID] += o[i];
		}
	}
	//statsContig.print(std::cout, CONTIG_ALL, this->getHeader().contigs_[ref_contigA].name, -1, ".");
	//statsContigPosition.print(std::cout, CONTIG_POS, this->getHeader().contigs_[ref_contigA].name, ref_positionA, ".");
	//statsContigPositionContig.print(std::cout, CONTIG_POS_CONTIG, this->getHeader().contigs_[ref_contigA].name, ref_positionA, this->getHeader().contigs_[ref_contigB].name);
	//statsContigSelf.print(std::cout, CONTIG_SELF, this->getHeader().contigs_[ref_contigA].name, -1, this->getHeader().contigs_[ref_contigA].name);

	contig_data.push_back(statsContig);
	contig_position_data.push_back(statsContigPosition);
	contig_position_contig_data.push_back(statsContigPositionContig);
	contig_self_data.push_back(statsContigSelf);

	/*
	for(U32 i = 0; i < this->getHeader().magic_.n_contigs; ++i){
		for(U32 j = i + 1; j < this->getHeader().magic_.n_contigs; ++j){
			if(pairwise_contigs[i][j].R.n_total == 0) continue;
			pairwise_contigs[i][j].print(std::cout, "contig-across", this->getHeader().contigs_[i].name, -1, this->getHeader().contigs_[j].name);
		}
	}
	*/
	std::cerr << contig_data.size() << ", " << contig_self_data.size() << ", " << contig_position_data.size() << ", " << contig_position_contig_data.size() << std::endl;

	return true;
}

bool TomahawkOutputReader::aggregate(support::aggregation_parameters& parameters){
	assert(this->index_ != nullptr);

	if(this->index_->isSorted() == false){
		std::cerr << "data is not sorted" << std::endl;
		return false;
	}

	if(this->index_->getMetaContainer().size() == 0){
		std::cerr << "Has no meta index..." << std::endl;
		return false;
	}

	typedef double (entry_type::*value_accessor_function)(void) const;
	typedef double (SummaryStatistics::*reduction_function)(void) const;
	typedef SummaryStatistics summary_statistics_type;

	value_accessor_function value_accessor = &entry_type::getR2;
	reduction_function reduction_accessor  = &SummaryStatistics::getMean;

	switch(parameters.aggregation_target){
	case(support::TWK_AGGREGATE_COUNT):   value_accessor = &entry_type::getR2; break;
	case(support::TWK_AGGREGATE_D):       value_accessor = &entry_type::getD;  break;
	case(support::TWK_AGGREGATE_DPrime):  value_accessor = &entry_type::getDPrime; break;
	case(support::TWK_AGGREGATE_R):       value_accessor = &entry_type::getR;  break;
	case(support::TWK_AGGREGATE_R2):      value_accessor = &entry_type::getR2; break;
	case(support::TWK_AGGREGATE_P1):      value_accessor = &entry_type::getP1; break;
	case(support::TWK_AGGREGATE_P2):      value_accessor = &entry_type::getP2; break;
	case(support::TWK_AGGREGATE_Q1):      value_accessor = &entry_type::getQ1; break;
	case(support::TWK_AGGREGATE_Q2):      value_accessor = &entry_type::getQ2; break;
	case(support::TWK_AGGREGATE_P_VALUE): value_accessor = &entry_type::getP;  break;
	case(support::TWK_AGGREGATE_LOG_P_VALUE): value_accessor = &entry_type::getLog10P;  break;
	}


	switch(parameters.reduction_target){
	case(support::TWK_AGGREGATE_REDUCE_COUNT): reduction_accessor = &summary_statistics_type::getCount; break;
	case(support::TWK_AGGREGATE_REDUCE_MEAN):  reduction_accessor = &summary_statistics_type::getMean;  break;
	case(support::TWK_AGGREGATE_REDUCE_MIN):   reduction_accessor = &summary_statistics_type::getMin;   break;
	case(support::TWK_AGGREGATE_REDUCE_MAX):   reduction_accessor = &summary_statistics_type::getMax;   break;
	case(support::TWK_AGGREGATE_REDUCE_SD):    reduction_accessor = &summary_statistics_type::getStandardDeviation; break;
	case(support::TWK_AGGREGATE_REDUCE_SUM):   reduction_accessor = &summary_statistics_type::getTotal; break;
	case(support::TWK_AGGREGATE_REDUCE_SUM_SQUARED): reduction_accessor = &summary_statistics_type::getTotalSquared; break;
	}

	// Projection
	// Setup matrix of summary statistics
	summary_statistics_type** matrix = new summary_statistics_type*[parameters.scene_x_pixels];
	for(U32 i = 0; i < parameters.scene_x_pixels; ++i)
		matrix[i] = new summary_statistics_type[parameters.scene_y_pixels];


	struct OffsetSupport{
	public:
		OffsetSupport() : range(0), min(0), max(0), cumulative(0){}
		OffsetSupport(const U64 range, const U64 min, const U64 max, const U64 cumulative):
			range(range),
			min(min),
			max(max),
			cumulative(cumulative)
		{

		}

		~OffsetSupport() = default;

	public:
		U64 range;
		U64 min;
		U64 max;
		U64 cumulative;
	};

	U64 cumulative_position = 0;
	OffsetSupport* cumulative_offsets = new OffsetSupport[this->getIndex().getMetaContainer().size()+1];
	cumulative_offsets[0]   = OffsetSupport(0,0,0,0);
	for(U32 i = 0; i < this->getHeader().magic_.n_contigs; ++i){
		U64 range = this->index_->getMetaContainer().at(i).max_position - this->index_->getMetaContainer().at(i).min_position;
		if(this->index_->getMetaContainer().at(i).max_position){
			++range;

			//cumulative_position += this->getHeader().contigs_[i].n_bases;
			//cumulative_position += this->index_->getMetaContainer().at(i).max_position + 1;
			cumulative_position += range;

			//std::cerr << range << ": " << this->index_->getMetaContainer().at(i).min_position << "->" << this->index_->getMetaContainer().at(i).max_position << std::endl;
		}

		cumulative_offsets[i].range = range;
		cumulative_offsets[i].min   = this->index_->getMetaContainer().at(i).min_position;
		cumulative_offsets[i].max   = this->index_->getMetaContainer().at(i).max_position;
		cumulative_offsets[i].range;
		cumulative_offsets[i+1].cumulative = cumulative_position;

		//std::cerr << "Cumulative: " << cumulative_position << std::endl;
	}

	//exit(1);

	std::cerr << helpers::timestamp("LOG") << "Pixel size: " << (U32)((double)cumulative_position/parameters.scene_x_pixels) << " bases..." << std::endl;

	// While there is blocks available
	while(this->nextBlock()){
		containers::OutputContainerReference o = this->getContainerReference();

		for(U32 i = 0; i < o.size(); ++i){
			//U32 fromBin = (U32)((double)(cumulative_offsets[o[i].AcontigID].cumulative + o[i].Aposition) / cumulative_position * parameters.scene_x_pixels);
			//U32 toBin   = (U32)((double)(cumulative_offsets[o[i].BcontigID].cumulative + o[i].Bposition) / cumulative_position * parameters.scene_y_pixels);

			U32 fromBin = (U32)((double)(cumulative_offsets[o[i].AcontigID].cumulative + o[i].Aposition - cumulative_offsets[o[i].AcontigID].min)  / cumulative_position * parameters.scene_x_pixels);
			U32 toBin   = (U32)((double)(cumulative_offsets[o[i].BcontigID].cumulative + o[i].Bposition - cumulative_offsets[o[i].BcontigID].min) / cumulative_position * parameters.scene_y_pixels);

			/*
			std::cerr << fromBin << "," << toBin << std::endl;
			std::cerr << cumulative_offsets[o[i].AcontigID].cumulative + o[i].Aposition - cumulative_offsets[o[i].AcontigID].min << std::endl;
			std::cerr << cumulative_offsets[o[i].BcontigID].cumulative + o[i].Bposition - cumulative_offsets[o[i].BcontigID].min << std::endl;
			std::cerr << "Min: " << cumulative_offsets[o[i].AcontigID].min << "," << cumulative_offsets[o[i].BcontigID].min << std::endl;
			std::cerr << "Cumsum: " << cumulative_offsets[o[i].AcontigID].cumulative << std::endl;
			*/

			// Target pixel : cumulative offset up to previous contig + current position / total cumulative offset


			//cumulative_offsets[o[i].AcontigID].min


			if(toBin >= parameters.scene_y_pixels){
				std::cerr << o[i].Aposition << "->" << cumulative_offsets[o[i].AcontigID].cumulative + o[i].Aposition << "/" << cumulative_position << "->" << (U64)((double)(cumulative_offsets[o[i].AcontigID].cumulative + o[i].Aposition) / cumulative_position * parameters.scene_x_pixels) << std::endl;
				std::cerr << o[i].Bposition << "->" << cumulative_offsets[o[i].BcontigID].cumulative + o[i].Bposition << "/" << cumulative_position << "->" << (U64)((double)(cumulative_offsets[o[i].BcontigID].cumulative + o[i].Bposition) / cumulative_position * parameters.scene_y_pixels) << std::endl;

				std::cerr << "Y: " << fromBin << "->" << toBin << std::endl;
				exit(1);
				toBin = parameters.scene_y_pixels - 1;
			}

			if(fromBin >= parameters.scene_x_pixels){
				std::cerr << o[i].Aposition << "->" << cumulative_offsets[o[i].AcontigID].cumulative + o[i].Aposition << "/" << cumulative_position << "->" << (U64)((double)(cumulative_offsets[o[i].AcontigID].cumulative + o[i].Aposition) / cumulative_position * parameters.scene_x_pixels) << std::endl;
				std::cerr << o[i].Bposition << "->" << cumulative_offsets[o[i].BcontigID].cumulative + o[i].Bposition << "/" << cumulative_position << "->" << (U64)((double)(cumulative_offsets[o[i].BcontigID].cumulative + o[i].Bposition) / cumulative_position * parameters.scene_y_pixels) << std::endl;

				std::cerr << "X: " << fromBin << "->" << toBin << std::endl;
				exit(1);
				fromBin = parameters.scene_x_pixels - 1;
			}


			//std::cerr << o[i].Aposition << "->" << cumulative_offsets[o[i].AcontigID] + o[i].Aposition - offset_min << "/" << cumulative_position << "->" << (U64)((double)(cumulative_offsets[o[i].AcontigID] + o[i].Aposition - offset_min) / cumulative_position * scene_x_length) << std::endl;
			//std::cerr << o[i].Bposition << "->" << cumulative_offsets[o[i].BcontigID] + o[i].Bposition - offset_min << "/" << cumulative_position << "->" << (U64)((double)(cumulative_offsets[o[i].BcontigID] + o[i].Bposition - offset_min) / cumulative_position * scene_y_length) << std::endl;

			//if(this->filters_.filter(o[i])) continue;
			//const double val = (o[i].*value_accessor)();
			matrix[fromBin][toBin] += (o[i].*value_accessor)();
			//std::cerr << val << "\t" << pow(val,2) << std::endl;
		}
	}

	for(U32 i = 0; i < parameters.scene_x_pixels; ++i){
		//std::cout << i << "/" << 0 << ":" << matrix[i][0].R2 << std::endl;
		matrix[i][0].calculate();
		// Mask values with a total number of observations <= 5
		//if(matrix[i][0].n_total > 1000)
		std::cout << (matrix[i][0].*reduction_accessor)();
		//else std::cout << 0;
		//std::cout << matrix[i][0].n_total;

		for(U32 j = 1; j < parameters.scene_y_pixels; ++j){
			matrix[i][j].calculate();
			// Mask values with a total number of observations <= 5
			//if(matrix[i][j].n_total > 1000)
			std::cout << '\t' << (matrix[i][j].*reduction_accessor)();
			//else std::cout << '\t' << 0;
			//std::cout << "\t" << matrix[i][j].n_total;
			//std::cout << i << "/" << j << ":" << matrix[i][j].R2 << std::endl;
		}
		std::cout.put('\n');
	}
	std::cout.flush();

	for(U32 i = 0; i < parameters.scene_x_pixels; ++i) delete [] matrix[i];
	delete [] matrix;
	delete [] cumulative_offsets;

	return true;
}

} /* namespace Tomahawk */
