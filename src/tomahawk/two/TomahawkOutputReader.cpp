#include "tomahawk/two/TomahawkOutputReader.h"

#include <bits/move.h>
#include <io/compression/gz_constants.h>
#include <io/compression/gz_header.h>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "support/helpers.h"
#include "tomahawk/two/TomahawkOutputStats.h"
#include "io/output_writer.h"

namespace tomahawk {

TomahawkOutputReader::TomahawkOutputReader() :
		filesize_(0),
		offset_end_of_data_(0),
		showHeader_(true),
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
	this->buffer_.deleteAll();
	this->data_.deleteAll();
	this->outputBuffer_.deleteAll();
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

	if(this->header_.magic_.major_version == 0 && this->header_.magic_.minor_version < 0.4){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Legacy file not supported..." << std::endl;
		return false;
	}

	this->offset_end_of_data_ = this->footer_.offset_end_of_data;

	return true;
}

int TomahawkOutputReader::parseBlock(const bool clear){
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

	// Reset compressed_buffer
	this->buffer_.reset();

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

	if(this->interval_tree_entries == nullptr)
		this->interval_tree_entries = new std::vector<interval_type>[this->getHeader().getMagic().getNumberContigs()];

	if(this->interval_tree == nullptr){
		this->interval_tree = new tree_type*[this->getHeader().getMagic().getNumberContigs()];
		for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i)
			this->interval_tree[i] = nullptr;
	}

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
	for(U32 i = 0; i < positions.size(); ++i){
		if(positions[i].find(',') != std::string::npos){
			std::vector<std::string> ret = helpers::split(positions[i], ',');
			if(ret.size() == 1){
				std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
				return false;

			} else if(ret.size() == 2){
				// parse left
				interval_type intervalLeft;
				if(this->__ParseRegionIndexed(ret[0], intervalLeft))
					this->interval_tree_entries[intervalLeft.contigID].push_back(interval_type(intervalLeft));

				// parse right
				interval_type intervalRight;
				if(this->__ParseRegionIndexed(ret[1], intervalRight))
					this->interval_tree_entries[intervalRight.contigID].push_back(interval_type(intervalRight));

			} else {
				std::cerr << helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
				return false;
			}
		}
		// Has no comma in string
		else {
			interval_type interval;
			if(this->__ParseRegionIndexed(positions[i], interval))
				this->interval_tree_entries[interval.contigID].push_back(interval_type(interval));
		}
	}

	return true;
}

bool TomahawkOutputReader::__ParseRegionIndexed(const std::string& region, interval_type& interval){
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
	//std::cerr << helpers::timestamp("LOG") << "Sorted: " << (int)this->getIndex().getController().isSorted << " partial: " << (int)this->getIndex().getController().isPartialSorted << std::endl;
	this->getHeader().getLiterals() += "\n##tomahawk_viewCommand=" + helpers::program_string();
	this->getHeader().getLiterals() += "\n##tomahawk_viewFilters=" + this->filters_.getInterpretedString() + " filter=NO regions=NO";

	//if(!this->OpenWriter())
	//	return false;

	if(this->showHeader_ == true){
		std::cout << this->getHeader().getLiterals() << '\n';
		std::cout << "FLAG\tCHROM_A\tPOS_A\tCHROM_B\tPOS_B\tREF_REF\tREF_ALT\tALT_REF\tALT_ALT\tD\tDprime\tR\tR2\tP\tChiSqModel\tChiSqTable\n";
	}

	// Natural output required parsing
	size_t n_total = 0;
	//if(this->writer_output_type == WRITER_TYPE::natural){
		while(this->parseBlock()){
			OutputContainerReference o = this->getContainerReference();
			n_total += o.size();
			for(U32 i = 0; i < o.size(); ++i){
				o[i].write(std::cout, this->getHeader().contigs_);
			}
		}
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
	this->getHeader().getLiterals() += "\n##tomahawk_viewCommand=" + helpers::program_string();
	if(this->filters_.any_filter_user_set)
		this->getHeader().getLiterals() += "\n##tomahawk_viewFilters=" + this->filters_.getInterpretedString() + " filter=YES regions=YES";

	if(this->showHeader_ == true){
		std::cout << this->getHeader().getLiterals() << '\n';
		std::cout << "FLAG\tCHROM_A\tPOS_A\tCHROM_B\tPOS_B\tREF_REF\tREF_ALT\tALT_REF\tALT_ALT\tD\tDprime\tR\tR2\tP\tChiSqModel\tChiSqTable\n";
	}

	if(this->interval_tree != nullptr){
		while(this->parseBlock()){
			output_container_reference_type o(this->data_);
			for(U32 i = 0; i < o.size(); ++i)
				this->__checkRegionNoIndex(o[i]);
		} // end while next block
	}

	return true;
}

bool TomahawkOutputReader::__viewFilter(void){
	this->getHeader().getLiterals() += "\n##tomahawk_viewCommand=" + helpers::program_string();
	this->getHeader().getLiterals() += "\n##tomahawk_viewFilters=" + this->filters_.getInterpretedString() + " filter=YES regions=NO";

	if(this->showHeader_ == true){
		std::cout << this->getHeader().getLiterals() << '\n';
		std::cout << "FLAG\tCHROM_A\tPOS_A\tCHROM_B\tPOS_B\tREF_REF\tREF_ALT\tALT_REF\tALT_ALT\tD\tDprime\tR\tR2\tP\tChiSqModel\tChiSqTable\n";
	}

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

bool TomahawkOutputReader::__checkRegionNoIndex(const entry_type& entry){
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
							entry.write(std::cout, this->getHeader().contigs_);
						}

						return true;
					} // end match
				} else { //  not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						entry.write(std::cout, this->getHeader().contigs_);
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
							entry.write(std::cout, this->getHeader().contigs_);
						}
						return true;
					} // end match
				} else { // not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						entry.write(std::cout, this->getHeader().contigs_);
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
		std::cerr << helpers::timestamp("LOG", "CONCAT") << "Opening input: " << files[0] << "..." << std::endl;

	if(!this->open(files[0])){
		std::cerr << helpers::timestamp("ERROR","TWO") << "Failed to parse: " << files[0] << "..." << std::endl;
		return false;
	}

	this->getHeader().getLiterals() += "\n##tomahawk_concatCommand=" + helpers::program_string();
	this->getHeader().getLiterals() += "\n##tomahawk_concatFiles=";
	for(U32 i = 0; i < files.size(); ++i)
		this->getHeader().getLiterals() += files[i] + ',';

	io::OutputWriter writer;
	if(!writer.open(output)){
		std::cerr << helpers::timestamp("ERROR","SORT") << "Failed to open: " << output << "..." << std::endl;
		return false;
	}
	writer.writeHeaders(this->getHeader());

	while(this->parseBlock())
		writer << this->data_;

	for(U32 i = 1; i < files.size(); ++i){
		if(!SILENT)
			std::cerr << helpers::timestamp("LOG", "CONCAT") << "Opening input: " << files[i] << "..." << std::endl;

		this->stream_.close();
		self_type second_reader;
		if(!second_reader.open(files[i])){
			std::cerr << helpers::timestamp("ERROR","TWO") << "Failed to parse: " << files[i] << "..." << std::endl;
			return false;
		}

		if(!(second_reader.getHeader() == this->getHeader())){
			std::cerr << "header mismatch" << std::endl;
		}

		while(second_reader.parseBlock())
			writer << this->data_;
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
		std::cerr << helpers::timestamp("ERROR","TWO") << "Failed to get file_list..." << std::endl;
		return false;
	}

	std::vector<std::string> files;
	std::string line;
	while(getline(file_list_read, line)){
		if(line.size() == 0){
			std::cerr << helpers::timestamp("WARNING","TWO") << "Empty line" << std::endl;
			break;
		}
		files.push_back(line);
	}

	return(this->__concat(files, output));
}

} /* namespace Tomahawk */
