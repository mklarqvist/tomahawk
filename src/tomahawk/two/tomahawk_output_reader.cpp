#include <io/output_writer.h>
#include <tomahawk/two/tomahawk_output_reader.h>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "io/compression/gz_constants.h"
#include "io/compression/gz_header.h"
#include "support/helpers.h"
#include "math/output_statistics.h"
#include "tomahawk_output_stats.h"

namespace tomahawk {

TomahawkOutputReader::TomahawkOutputReader() :
		filesize_(0),
		offset_end_of_data_(0),
		showHeader_(true),
		output_json_(false),
		index_(nullptr),
		buffer_(3000000),
		data_(3000000),
		outputBuffer_(3000000),
		interval_tree(nullptr),
		interval_tree_entries(nullptr)
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
}

bool TomahawkOutputReader::open(const std::string input){
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
	if(this->showHeader_ == true){
		std::cout << this->getHeader().getLiterals() << '\n';
		std::cout << "FLAG\tCHROM_A\tPOS_A\tCHROM_B\tPOS_B\tREF_REF\tREF_ALT\tALT_REF\tALT_ALT\tD\tDprime\tR\tR2\tP\tChiSqModel\tChiSqTable\n";
	}

	return true;
}

int TomahawkOutputReader::parseBlock(const bool clear, const bool clear_raw){
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

OutputContainer TomahawkOutputReader::getContainerVariants(const U64 n_variants){
	size_t n_variants_loaded = 0;
	this->data_.reset();
	this->data_.resize(n_variants*sizeof(entry_type) + 65536); // make room for data
	while(true){
		if(!this->parseBlock(false))
			break;

		n_variants_loaded = this->data_.size() / sizeof(entry_type);
		//std::cerr << n_variants_loaded << "/" << n_variants << '\t' << this->data_.size() << std::endl;
		if(n_variants_loaded >= n_variants)
			break;
	}

	return(OutputContainer(this->data_));
}

OutputContainer TomahawkOutputReader::getContainerBytes(const size_t l_data){
	const U64 start_position = this->stream_.tellg();
	this->data_.reset();
	this->data_.resize(l_data + 65536); // make room for data
	U64 data_loaded = 0;
	while(true){
		if(!this->parseBlock(false))
			break;

		data_loaded = (U64)this->stream_.tellg() - start_position;
		if(data_loaded >= l_data)
			break;

	}

	return(OutputContainer(this->data_));
}

OutputContainer TomahawkOutputReader::getContainerBlocks(const U32 n_blocks){
	this->data_.reset();
	for(U32 i = 0; i < n_blocks; ++i){
		if(!this->parseBlock(false))
			break;
	}

	return(OutputContainer(this->data_));
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


OutputContainerReference TomahawkOutputReader::getContainerReferenceBlock(const U32 blockID){
	if(!this->seekBlock(blockID)){
		this->parseBlock();
		return(OutputContainerReference());
	} else {
		this->parseBlock();
		return(OutputContainerReference(this->data_));
	}
}

OutputContainer TomahawkOutputReader::getContainerBlock(const U32 blockID){
	if(!this->seekBlock(blockID)){
		this->parseBlock();
		return(OutputContainer());
	} else {
		this->parseBlock();
		return(OutputContainer(this->data_));
	}
}

OutputContainerReference TomahawkOutputReader::getContainerReferenceBlock(std::vector<U32> blocks){
	if(!this->seekBlock(blocks[0])){
		this->parseBlock();
		return(OutputContainerReference());
	} else {
		for(U32 i = 0; i < blocks.size(); ++i){
			if(!this->parseBlock(false))
				break;
		}
		return(OutputContainerReference(this->data_));
	}
}

OutputContainer TomahawkOutputReader::getContainerBlock(std::vector<U32> blocks){
	if(!this->seekBlock(blocks[0])){
		this->parseBlock();
		return(OutputContainer());
	} else {
		for(U32 i = 0; i < blocks.size(); ++i){
			if(!this->parseBlock(false))
				break;
		}
		return(OutputContainer(this->data_));
	}
}

bool TomahawkOutputReader::addRegions(std::vector<std::string>& positions){
	if(positions.size() == 0)
		return true;

	// No contigs (usually means data has not been loaded)
	if(this->getHeader().getMagic().getNumberContigs() == 0)
		return false;

	// Construct interval tree and interval vector if not set
	if(this->interval_tree_entries == nullptr)
		this->interval_tree_entries = new std::vector<interval_type>[this->getHeader().getMagic().getNumberContigs()];

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
				if(this->__ParseRegion(ret[0], intervalLeft))
					this->interval_tree_entries[intervalLeft.contigID].push_back(interval_type(intervalLeft));

				// parse right
				interval_type intervalRight;
				if(this->__ParseRegion(ret[1], intervalRight))
					this->interval_tree_entries[intervalRight.contigID].push_back(interval_type(intervalRight));

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

	if(this->showHeader_ == true)
		this->printHeader(std::cout, extra);

	const std::string version_string = std::to_string(this->header_.magic_.major_version) + "." + std::to_string(this->header_.magic_.minor_version) + "." + std::to_string(this->header_.magic_.patch_version);


	// Natural output required parsing
	size_t n_total = 0;
	//if(this->writer_output_type == WRITER_TYPE::natural){
	typedef std::ostream& (entry_type::*func)(std::ostream& os, const contig_type* const contigs) const;
	func a = &entry_type::write;
	if(this->output_json_){
		a = &entry_type::writeJSON;
		std::cout << "{\n\"type\":\"tomahawk\",\n\"sorted\":" << (this->getIndex().getController().isSorted ? "true" : "false") << ",\n\"partial_sort\":" << (this->getIndex().getController().isPartialSorted ? "true" : "false");
		std::cout <<  ",\n\"version\":\"" << version_string << "\",\n\"data\":[\n";
	}

	//this->parseBlock();
		while(this->parseBlock()){
			OutputContainerReference o = this->getContainerReference();
			n_total += o.size();
			//if(o.size() == 0) break;

			(o[0].*a)(std::cout, this->getHeader().contigs_);
			for(U32 i = 1; i < o.size(); ++i){
				//std::cout << ",\n";
				(o[i].*a)(std::cout, this->getHeader().contigs_);
			}
		}

		//std::cout << "]\n}\n";
		//std::cerr << "total: " << n_total << std::endl;
	//}
	// Binary output without filtering simply writes it back out
/*
	else if(this->writer_output_type == WRITER_TYPE::binary){
		while(this->parseBlock()){
			OutputContainerReference o(this->compressed_buffer);
			//this->writer->write(this->data_);
			std::cout << o[0] << std::endl;
		}
	}
*/
	return true;
}

bool TomahawkOutputReader::__viewRegion(void){
	std::vector<std::string> extra;
	extra.push_back("##tomahawk_viewCommand=" + helpers::program_string());
	if(this->filters_.any_filter_user_set)
		extra.push_back("##tomahawk_viewFilters=" + this->filters_.getInterpretedString() + " filter=YES regions=YES");

	if(this->showHeader_ == true)
		this->printHeader(std::cout, extra);

	if(this->interval_tree == nullptr)
		return false;

	typedef bool (TomahawkOutputReader::*func_slice)(const entry_type& entry);
	func_slice region_function = &TomahawkOutputReader::__checkRegionUnsorted;
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

		while(this->parseBlock()){
			output_container_reference_type o(this->data_);
			for(U32 i = 0; i < o.size(); ++i)
				(this->*region_function)(o[i]);
		}
	}

	return true;
}

bool TomahawkOutputReader::__viewFilter(void){
	std::vector<std::string> extra;
	extra.push_back("##tomahawk_viewCommand=" + helpers::program_string());
	extra.push_back("##tomahawk_viewFilters=" + this->filters_.getInterpretedString() + " filter=YES regions=NO");

	if(this->showHeader_ == true)
		this->printHeader(std::cout, extra);

	while(this->parseBlock()){
		output_container_reference_type o(this->data_);
		for(U32 i = 0; i < o.size(); ++i){
			if(this->filters_.filter(o[i])){
				o[i].write(std::cout, this->getHeader().contigs_);
			}
		}
	} // end while next block
	return true;
}

bool TomahawkOutputReader::__checkRegionSorted(const entry_type& entry){
	typedef std::ostream& (entry_type::*func)(std::ostream& os, const contig_type* const contigs) const;
	func a = &entry_type::write;

	// If iTree for contigA exists
	if(this->interval_tree[entry.AcontigID] != nullptr){
		std::vector<interval_type> rets = this->interval_tree[entry.AcontigID]->findOverlapping(entry.Aposition, entry.Aposition);
		if(rets.size() > 0){
			for(U32 i = 0; i < rets.size(); ++i){
				if(rets[i].value != nullptr){ // if linked
					if((entry.BcontigID == rets[i].value->contigID) &&
					   (entry.Bposition >= rets[i].value->start &&
					    entry.Bposition <= rets[i].value->stop)){
						if(this->filters_.filter(entry)){
							//entry.write(std::cout, this->contigs);
							//*this->writer << entry;
							(entry.*a)(std::cout, this->getHeader().contigs_);
						}

						return true;
					} // end match
				} else { //  not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						(entry.*a)(std::cout, this->getHeader().contigs_);
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
	func a = &entry_type::write;

	// If iTree for contigA exists
	if(this->interval_tree[entry.AcontigID] != nullptr){
		std::vector<interval_type> rets = this->interval_tree[entry.AcontigID]->findOverlapping(entry.Aposition, entry.Aposition);
		if(rets.size() > 0){
			for(U32 i = 0; i < rets.size(); ++i){
				if(rets[i].value != nullptr){ // if linked
					if((entry.BcontigID == rets[i].value->contigID) &&
					   (entry.Bposition >= rets[i].value->start &&
					    entry.Bposition <= rets[i].value->stop)){
						if(this->filters_.filter(entry)){
							//entry.write(std::cout, this->contigs);
							//*this->writer << entry;
							(entry.*a)(std::cout, this->getHeader().contigs_);
						}

						return true;
					} // end match
				} else { //  not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						(entry.*a)(std::cout, this->getHeader().contigs_);
					}

					return true;
				}
			}
		}
	}

	// If iTree for contigB exists
	if(this->interval_tree[entry.BcontigID] != nullptr){
		std::vector<interval_type> rets = this->interval_tree[entry.BcontigID]->findOverlapping(entry.Bposition, entry.Bposition);
		if(rets.size() > 0){
			for(U32 i = 0; i < rets.size(); ++i){
				if(rets[i].value != nullptr){ // if linked
					if((entry.AcontigID == rets[i].value->contigID) &&
					   (entry.Aposition >= rets[i].value->start &&
						entry.Aposition <= rets[i].value->stop)){
						if(this->filters_.filter(entry)){
							//entry.write(std::cout, this->contigs);
							//*this->writer << entry;
							(entry.*a)(std::cout, this->getHeader().contigs_);
						}
						return true;
					} // end match
				} else { // not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						(entry.*a)(std::cout, this->getHeader().contigs_);
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

	io::OutputWriterFile writer;
	if(!writer.open(output)){
		std::cerr << helpers::timestamp("ERROR","SORT") << "Failed to open: " << output << "..." << std::endl;
		return false;
	}
	writer.writeHeaders(this->getHeader());

	while(this->parseBlock()){
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

		while(second_reader.parseBlock()){
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

	if(this->parseBlock() == false){
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
	OutputContainerReference o = this->getContainerReference();
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
	while(this->parseBlock()){
		OutputContainerReference o = this->getContainerReference();

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

} /* namespace Tomahawk */
