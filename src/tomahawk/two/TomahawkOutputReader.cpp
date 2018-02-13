#include "../two/TomahawkOutputReader.h"

#include <bits/move.h>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "../../io/compression/GZFConstants.h"
#include "../../io/compression/GZFHeader.h"
#include "../../support/helpers.h"
#include "../two/TomahawkOutputStats.h"

namespace Tomahawk {

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
		for(U32 i = 0; i < this->getHeader().magic_.getNumberContigs(); ++i){
			delete this->interval_tree[i];
		}
		delete [] this->interval_tree;
	}
}

bool TomahawkOutputReader::open(const std::string input){
	if(input.size() == 0){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "No input filename..." << std::endl;
		return false;
	}

	this->stream_.open(input, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to open file handle: " << input << std::endl;
	}
	this->filesize_ = this->stream_.tellg();

	this->stream_.seekg(this->filesize_ - TWK_FOOTER_LENGTH);
	this->stream_ >> this->footer_;
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Stream corrupted after loading footer..." << std::endl;
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

	return true;
}

int TomahawkOutputReader::parseBlock(const bool clear){
	// Stream died
	if(!this->stream_.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Stream died!" << std::endl;
		return -1;
	}

	// EOF
	// tellg will always return a positive value here
	// or it would've failed at good() check
	if((U64)this->stream_.tellg() == this->offset_end_of_data_)
		return 0;

	// Read TGZF header
	this->buffer_.resize(sizeof(tgzf_header_type));
	this->stream_.read(this->buffer_.data(),  IO::Constants::TGZF_BLOCK_HEADER_LENGTH);
	const tgzf_header_type* h = reinterpret_cast<const tgzf_header_type*>(this->buffer_.data());
	this->buffer_.n_chars = IO::Constants::TGZF_BLOCK_HEADER_LENGTH;
	if(!h->Validate()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to validate!" << std::endl;
		return -2;
	}

	this->buffer_.resize(h->BSIZE); // make sure all data will fit

	// Recast because if compressed_buffer is (actually) resized then the pointer address is incorrect
	// resulting in segfault
	h = reinterpret_cast<const tgzf_header_type*>(this->buffer_.data());

	this->stream_.read(&this->buffer_.buffer[IO::Constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - IO::Constants::TGZF_BLOCK_HEADER_LENGTH);
	if(!this->stream_.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Truncated file..." << std::endl;
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
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed inflate!" << std::endl;
		return -4;
	}

	if(this->data_.size() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Empty data!" << std::endl;
		return 0;
	}

	// Reset compressed_buffer
	this->buffer_.reset();

	// Reset iterator position and size
	//this->iterator_position_block = 0;
	//this->size = this->data_.size() / sizeof(entry_type);

	// Validity check
	if(this->data_.size() % sizeof(entry_type) != 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Data is corrupted!" << std::endl;
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

bool TomahawkOutputReader::seekBlock(const U32 blockID){
	if(blockID > this->getIndex().getContainer().size()){
		std::cerr << Helpers::timestamp("ERROR","TOI") << "Illegal blockID (" << blockID << ">" << this->getIndex().getContainer().size() << ")!" << std::endl;
		return false;
	}

	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Stream is bad!" << std::endl;
		return false;
	}

	this->stream_.seekg(this->getIndex().getContainer()[blockID].byte_offset);
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Stream is bad following seek!" << std::endl;
		return false;
	}

	return(this->parseBlock());
}

OutputContainerReference TomahawkOutputReader::getContainerReferenceBlock(const U32 blockID){
	if(!this->seekBlock(blockID)){
		return(OutputContainerReference());
	} else {
		this->parseBlock();
		return(OutputContainerReference(this->data_));
	}
}

OutputContainer TomahawkOutputReader::getContainerBlock(const U32 blockID){
	if(!this->seekBlock(blockID)){
		return(OutputContainer());
	} else {
		this->parseBlock();
		return(OutputContainer(this->data_));
	}
}

OutputContainerReference TomahawkOutputReader::getContainerReferenceBlock(std::vector<U32> blocks){
	if(!this->seekBlock(blocks[0])){
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
		this->interval_tree_entries = new std::vector<interval_type>[this->getHeader().magic_.getNumberContigs()];

	if(this->interval_tree == nullptr){
		this->interval_tree = new tree_type*[this->getHeader().magic_.getNumberContigs()];
		for(U32 i = 0; i < this->getHeader().magic_.getNumberContigs(); ++i)
			this->interval_tree[i] = nullptr;
	}

	if(!this->__addRegions(positions))
		return false;

	for(U32 i = 0; i < this->getHeader().magic_.getNumberContigs(); ++i){
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
			std::vector<std::string> ret = Helpers::split(positions[i], ',');
			if(ret.size() == 1){
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
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
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
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
	std::vector<std::string> ret = Helpers::split(region, ':');

	// If vector does not contain a colon
	if(ret.size() == 1){
		if(ret[0].find('-') != std::string::npos){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}

		// is contigID only
		const S32 contigID = this->getHeader().getContigID(ret[0]);
		if(contigID < 0){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << region << " is not defined in the header!" << std::endl;
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
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << ret[0] << " is not defined in the header!" << std::endl;
			return false;
		}

		std::vector<std::string> retPos = Helpers::split(ret[1], '-');
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
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}
	}
	// contains > 1 colons
	// illegal
	else {
		std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::view(void){
	if(this->interval_tree != nullptr) // If regions have been set: use region-filter function
		return(this->__viewRegion());
	//else if(this->filter.any_filter_user_set){
	//	return(this->__viewFilter()); // Otherwise normal filter function
	 else
		return(this->__viewOnly());
}

bool TomahawkOutputReader::__viewOnly(void){
	std::cerr << Helpers::timestamp("LOG") << "Sorted: " << (int)this->getIndex().getController().isSorted << " partial: " << (int)this->getIndex().getController().isPartialSorted << std::endl;
	this->getHeader().getLiterals() += "\n##tomahawk_viewCommand=" + Helpers::program_string();
	this->getHeader().getLiterals() += "\n##tomahawk_viewFilters=" + this->filters_.getInterpretedString() + " filter=NO regions=FALSE";

	//if(!this->OpenWriter())
	//	return false;

	if(this->showHeader_ == true){
		std::cout << this->getHeader().getLiterals() << '\n';
	}

	// Natural output required parsing
	size_t n_total = 0;
	//if(this->writer_output_type == WRITER_TYPE::natural){
		while(this->parseBlock()){
			OutputContainerReference o = this->getContainerReference();
			n_total += o.size();
			for(U32 i = 0; i < o.size(); ++i)
				std::cout << o[i] << '\n';
		}
		std::cerr << "total: " << n_total << std::endl;
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
	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "TWO") << "Unindexed query..." << std::endl;

	if(this->interval_tree != nullptr){
		const entry_type*  entry = nullptr;


		while(this->parseBlock()){
			output_container_reference_type o(this->data_);
			for(U32 i = 0; i < o.size(); ++i)
				this->__checkRegionNoIndex(o[i]);
		} // end while next variant
	}

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
					   (entry.Bposition >= rets[i].value->start && entry.Bposition <= rets[i].value->stop)){
						if(this->filters_.filter(entry)){
							//entry.write(std::cout, this->contigs);
							//*this->writer << entry;
							std::cout << entry << '\n';
						}

						return true;
					} // end match
				} else { //  not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						std::cout << entry << '\n';
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
					   (entry.Aposition >= rets[i].value->start && entry.Aposition <= rets[i].value->stop)){
						if(this->filters_.filter(entry)){
							//entry.write(std::cout, this->contigs);
							//*this->writer << entry;
							std::cout << entry << '\n';
						}
						return true;
					} // end match
				} else { // not linked
					if(this->filters_.filter(entry)){
						//entry.write(std::cout, this->contigs);
						//*this->writer << entry;
						std::cout << entry << '\n';
					}

					return true;
				}
			}
		} // end if any hit in iTree b
	} // end iTree b

	return false;
}

/*
bool TomahawkOutputReader::OpenWriter(void){
	if(this->writer_output_type == WRITER_TYPE::natural){
		this->writer = new OutputWriterNatural(this->contigs, &this->header);
	}
	else this->writer = new OutputWriter(this->contigs, &this->header);

	if(!this->writer->open())
		return false;

	if(this->output_header)
		this->writer->writeHeader(this->literals);

	return true;
}

bool TomahawkOutputReader::OpenWriter(const std::string output_file){
	if(this->writer_output_type == WRITER_TYPE::natural){
		this->writer = new OutputWriterNatural(this->contigs, &this->header);
	}
	else this->writer = new OutputWriter(this->contigs, &this->header);

	if(!this->writer->open(output_file)){
		std::cerr << Helpers::timestamp("ERROR","WRITER") << "Failed to open output file: " << output_file << std::endl;
		return false;
	}

	if(this->output_header)
		this->writer->writeHeader(this->literals);

	return true;
}

bool TomahawkOutputReader::__viewRegion(void){
	if(!this->OpenWriter())
		return false;

	// If indexed and expanded
	if(this->toi_reader.ERROR_STATE == toi_reader_type::TOI_OK && (this->toi_reader.getIsSortedExpanded())){
		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "TWO") << "Indexed query..." << std::endl;
		return(this->__viewRegionIndexed());
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "TWO") << "Unindexed query..." << std::endl;

	if(this->interval_tree != nullptr){
		const entry_type*  entry = nullptr;


		while(this->parseBlock()){
			output_container_reference_type o(this->data_);
			for(U32 i = 0; i < o.size(); ++i)
				this->__checkRegionNoIndex(o[i]);
		} // end while next variant
	}

	return true;
}
*/

bool TomahawkOutputReader::__viewRegionIndexed(void){
	/*
	if(this->interval_tree == nullptr){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Interval tree not set!" << std::endl;
		return false;
	}

	if(!this->__ParseRegionIndexedBlocks()){
		std::cerr << Helpers::timestamp("LOG","TOI") << "No valid entries..." << std::endl;
		return false;
	}

	// Init
	const entry_type* two_entry = nullptr;

	// Todo
	// sort entries
	// merge
	// for i in entries: seek, uncompress, and jump or limit
	for(U32 i = 0; i < this->interval_totempole_enties->size(); ++i){
		const totempole_sorted_entry_type& entry = this->interval_totempole_enties->at(i);
		const U32 block_length = entry.toBlock - entry.fromBlock;

		// 1 entry
		if(block_length == 0){
			if(!this->seekBlock(entry.fromBlock)){
				std::cerr << Helpers::timestamp("ERROR","TWO") << "Could not get block" << std::endl;
				return false;
			}
			this->iterator_position_block = entry.fromBlock_entries_offset;


			while(this->nextVariantLimited(two_entry)){
				this->__checkRegionIndex(two_entry);
			} // end while next variant
		}
		// 2 entries
		else if(block_length == 1){
			// First one
			if(!this->seekBlock(entry.fromBlock)){
				std::cerr << Helpers::timestamp("ERROR","TWO") << "Could not get block" << std::endl;
				return false;
			}
			this->iterator_position_block = entry.fromBlock_entries_offset;

			while(this->nextVariantLimited(two_entry)){
				this->__checkRegionIndex(two_entry);
			} // end while next variant

			// Second one
			if(!this->seekBlock(entry.toBlock)){
				std::cerr << Helpers::timestamp("ERROR","TWO") << "Could not get block" << std::endl;
				return false;
			}
			//this->position = entry.toBlock_entries_offset;

			while(this->nextVariantLimited(two_entry)){
				this->__checkRegionIndex(two_entry);
			} // end while next variant

		}
		// > 2 entries
		else {
			// First block
			U32 j = entry.fromBlock;
			if(!this->seekBlock(j)){
				std::cerr << Helpers::timestamp("ERROR","TWO") << "Could not get block" << std::endl;
				return false;
			}
			this->iterator_position_block = entry.fromBlock_entries_offset;

			while(this->nextVariantLimited(two_entry)){
				this->__checkRegionIndex(two_entry);
			} // end while next variant
			++j;

			// Middle blocks
			for(; j < entry.toBlock - 1; ++j){
				if(!this->seekBlock(j)){
					std::cerr << Helpers::timestamp("ERROR","TWO") << "Could not get block" << std::endl;
					return false;
				}

				while(this->nextVariantLimited(two_entry)){
					this->__checkRegionIndex(two_entry);
				} // end while next variant
			}

			// last block
			if(!this->seekBlock(j)){
				std::cerr << Helpers::timestamp("ERROR","TWO") << "Could not get block" << std::endl;
				return false;
			}
			//this->position = entry.toBlock_entries_offset;

			while(this->nextVariantLimited(two_entry)){
				this->__checkRegionIndex(two_entry);
			} // end while next variant
		}
	}

	 */
	return true;
}

/*
bool TomahawkOutputReader::__checkRegionIndex(const entry_type& entry){
	// If iTree for contigA exists
	if(this->interval_tree[entry.AcontigID] != nullptr){
		std::vector<interval_type> rets = this->interval_tree[entry.AcontigID]->findOverlapping(entry.Aposition, entry.Aposition);
		if(rets.size() > 0){
			for(U32 i = 0; i < rets.size(); ++i){
				if(this->filter.filter(entry))
					*this->writer << entry;

				return true;
			}
		}
	}
	return false;
}

bool TomahawkOutputReader::__checkRegionNoIndex(const entry_type& entry){
	// If iTree for contigA exists
	if(this->interval_tree[entry.AcontigID] != nullptr){
		std::vector<interval_type> rets = this->interval_tree[entry.AcontigID]->findOverlapping(entry.Aposition, entry.Aposition);
		if(rets.size() > 0){
			for(U32 i = 0; i < rets.size(); ++i){
				if(rets[i].value != nullptr){ // if linked
					if((entry.BcontigID == rets[i].value->contigID) &&
					   (entry.Bposition >= rets[i].value->start && entry.Bposition <= rets[i].value->stop)){
						if(this->filter.filter(entry))
							//entry.write(std::cout, this->contigs);
							*this->writer << entry;

						return true;
					} // end match
				} else { //  not linked
					if(this->filter.filter(entry))
						//entry.write(std::cout, this->contigs);
						*this->writer << entry;

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
					   (entry.Aposition >= rets[i].value->start && entry.Aposition <= rets[i].value->stop)){
						if(this->filter.filter(entry)){
							//entry.write(std::cout, this->contigs);
							*this->writer << entry;
						}
						return true;
					} // end match
				} else { // not linked
					if(this->filter.filter(entry))
						//entry.write(std::cout, this->contigs);
						*this->writer << entry;

					return true;
				}
			}
		} // end if any hit in iTree b
	} // end iTree b

	return false;
}

bool TomahawkOutputReader::__viewOnly(void){
	this->literals += "\n##tomahawk_viewCommand=" + Helpers::program_string();
	this->literals += "\n##tomahawk_viewFilters=" + this->filter.getInterpretedString() + " filter=NO regions=FALSE";

	if(!this->OpenWriter())
		return false;

	// Natural output required parsing
	size_t n_total = 0;
	if(this->writer_output_type == WRITER_TYPE::natural){
		while(this->parseBlock()){
			OutputContainerReference o = this->getContainerReference();
			std::cerr << o.size() << '\t' << this->data_.size() << std::endl;
			n_total += o.size();
			for(U32 i = 0; i < o.size(); ++i)
				std::cout << o[i] << '\n';
		}
		std::cerr << "total: " << n_total << std::endl;
	}
	// Binary output without filtering simply writes it back out
	else if(this->writer_output_type == WRITER_TYPE::binary){
		while(this->parseBlock()){
			OutputContainerReference o(this->compressed_buffer);
			//this->writer->write(this->data_);
			std::cout << o[0] << std::endl;
		}
	}

	return true;
}

bool TomahawkOutputReader::__viewFilter(void){
	this->literals += "\n##tomahawk_viewCommand=" + Helpers::program_string();
	this->literals += "\n##tomahawk_viewFilters=" + this->filter.getInterpretedString() + " filter=YES regions=FALSE";

	if(!this->OpenWriter())
		return false;

	while(this->parseBlock()){
		output_container_reference_type o(this->data_);
		for(U32 i = 0; i < o.size(); ++i){
			if(this->filter.filter(o[i]))
				*this->writer << o[i];
		}
	} // end while next variant
	return true;
}

bool TomahawkOutputReader::addRegionsIndexed(std::vector<std::string>& positions){
	for(U32 i = 0; i < positions.size(); ++i){
		if(positions[i].find(',') != std::string::npos){
			std::vector<std::string> ret = Helpers::split(positions[i], ',');
			if(ret.size() == 1){
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
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
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
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

bool TomahawkOutputReader::addRegionsUnindexed(std::vector<std::string>& positions){
	for(U32 i = 0; i < positions.size(); ++i){
		// Pattern cA:pAf-pAt;cB:pBf-pBt
		if(positions[i].find(',') != std::string::npos){
			//std::cerr << "linked intervals" << std::endl;
			std::vector<std::string> ret = Helpers::split(positions[i], ',');
			if(ret.size() == 1){
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
				return false;

			} else if(ret.size() == 2){
				// parse left
				interval_type intervalLeft;
				if(!this->__ParseRegion(ret[0], intervalLeft))
					return false;

				// parse right
				interval_type intervalRight;
				if(!this->__ParseRegion(ret[1], intervalRight))
					return false;

				// Todo: WARNING
				// This results in illegal pointers if the vector resizes
				// and pointers change
				this->interval_tree_entries[intervalLeft.contigID].push_back(interval_type(intervalLeft));
				this->interval_tree_entries[intervalRight.contigID].push_back(interval_type(intervalRight));
				if(intervalLeft.contigID != intervalRight.contigID){
					this->interval_tree_entries[intervalLeft.contigID].back().value  = &this->interval_tree_entries[intervalRight.contigID].back();
					this->interval_tree_entries[intervalRight.contigID].back().value = &this->interval_tree_entries[intervalLeft.contigID].back();
				} else {
					this->interval_tree_entries[intervalLeft.contigID].back().value = &this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2];
					this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2].value = &this->interval_tree_entries[intervalLeft.contigID].back();
				}

			} else {
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
				return false;
			}
		} else {
			interval_type interval;
			if(!this->__ParseRegion(positions[i], interval))
				return false;

			this->interval_tree_entries[interval.contigID].push_back(interval_type(interval));
		}
	}
	return true;
}

bool TomahawkOutputReader::addRegions(std::vector<std::string>& positions){
	if(positions.size() == 0)
		return true;

	if(this->interval_tree_entries == nullptr)
		this->interval_tree_entries = new std::vector<interval_type>[this->header.n_contig];

	if(this->interval_tree == nullptr){
		this->interval_tree = new tree_type*[this->header.n_contig];
		for(U32 i = 0; i < this->header.n_contig; ++i)
			this->interval_tree[i] = nullptr;
	}

	if(this->toi_reader.ERROR_STATE == toi_reader_type::TOI_OK && (this->toi_reader.getIsSortedExpanded())){
		if(!this->addRegionsIndexed(positions))
			return false;
	} else {
		if(!this->addRegionsUnindexed(positions))
			return false;
	}

	for(U32 i = 0; i < this->header.n_contig; ++i){
		if(this->interval_tree_entries[i].size() != 0){
			this->interval_tree[i] = new tree_type(this->interval_tree_entries[i]);
		} else
			this->interval_tree[i] = nullptr;
	}

	return true;
}

bool TomahawkOutputReader::__ParseRegion(const std::string& region, interval_type& interval){
	std::vector<std::string> ret = Helpers::split(region, ':');
	if(ret.size() == 1){
		if(ret[0].find('-') != std::string::npos){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}

		// is contigID only
		U32* contigID;
		if(!this->contig_htable->GetItem(&region[0], &region, contigID, region.size())){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << region << " is not defined in the header!" << std::endl;
			return false;
		}
		interval(*contigID, 0, this->contigs[*contigID].n_bases);

	} else if(ret.size() == 2){
		// is contigID:pos-pos
		U32* contigID;
		if(!this->contig_htable->GetItem(&ret[0][0], &ret[0], contigID, ret[0].size())){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << ret[0] << " is not defined in the header!" << std::endl;
			return false;
		}

		std::vector<std::string> retPos = Helpers::split(ret[1], '-');
		if(retPos.size() == 1){
			// only one pos
			const double pos = std::stod(retPos[0]);
			interval(*contigID, pos, pos);

		} else if(retPos.size() == 2){
			// is two positions
			double posA = std::stod(retPos[0]);
			double posB = std::stod(retPos[1]);

			if(posB < posA)
				std::swap(posA, posB);

			interval(*contigID, posA, posB);

		} else {
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}
	} else {
		std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::__ParseRegionIndexed(const std::string& region, interval_type& interval){
	std::vector<std::string> ret = Helpers::split(region, ':');

	// If vector does not contain a colon
	if(ret.size() == 1){
		if(ret[0].find('-') != std::string::npos){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}

		// is contigID only
		U32* contigID;
		if(!this->contig_htable->GetItem(&region[0], &region, contigID, region.size())){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << region << " is not defined in the header!" << std::endl;
			return false;
		}
		interval(*contigID, 0, this->contigs[*contigID].n_bases);
		interval.state = interval_type::INTERVAL_TYPE::INTERVAL_CONTIG_ONLY;
	}
	// If vector contain colon
	else if(ret.size() == 2){
		// is contigID:pos-pos
		U32* contigID;
		if(!this->contig_htable->GetItem(&ret[0][0], &ret[0], contigID, ret[0].size())){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << ret[0] << " is not defined in the header!" << std::endl;
			return false;
		}

		std::vector<std::string> retPos = Helpers::split(ret[1], '-');
		if(retPos.size() == 1){
			// only one pos
			const double pos = std::stod(retPos[0]);
			//std::cerr << "single position: " << pos << std::endl;
			interval(*contigID, pos, pos);
			interval.state = interval_type::INTERVAL_TYPE::INTERVAL_POSITION;

		} else if(retPos.size() == 2){
			// is two positions
			double posA = std::stod(retPos[0]);
			double posB = std::stod(retPos[1]);

			// Swap pA and pB iff pB > pA
			if(posB < posA)
				std::swap(posA, posB);

			interval(*contigID, posA, posB);
			interval.state = interval_type::INTERVAL_TYPE::INTERVAL_FULL;

		} else {
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}
	}
	// contains > 1 colons
	// illegal
	else {
		std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::__ParseRegionIndexedBlocks(void){
	if(this->interval_tree_entries == nullptr){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "No data is set" << std::endl;
		return false;
	}

	if(this->interval_totempole_enties == nullptr)
		this->interval_totempole_enties = new std::vector<totempole_sorted_entry_type>;

	for(U32 k = 0; k < this->header.n_contig; ++k){
		for(U32 i = 0; i < this->interval_tree_entries[k].size(); ++i){
			const interval_type& interval = this->interval_tree_entries[k][i];
			if(interval.state == interval_type::INTERVAL_TYPE::INTERVAL_CONTIG_ONLY){
				// Contig only
				//std::cerr << "contig only: " << interval << std::endl;
				totempole_sorted_entry_type entry;
				if(!this->toi_reader.findOverlap(interval.contigID, entry)){
					//std::cerr << "could not find: " << interval << std::endl;
					continue;
				}
				//std::cerr << "contigID found: " << entry << std::endl;
				this->interval_totempole_enties->push_back(entry);

			} else if(interval.state == interval_type::INTERVAL_TYPE::INTERVAL_POSITION){
				//std::cerr << "contig:posiiton only: " << interval << std::endl;
				totempole_sorted_entry_type entry;
				if(!this->toi_reader.findOverlap(interval.contigID, interval.start, entry)){
					//std::cerr << "could not find: " << interval << std::endl;
					continue;
				}
				//std::cerr << "contigID:pos found: " << entry << std::endl;
				this->interval_totempole_enties->push_back(entry);

			} else {
				//std::cerr << "full interval: " << interval << std::endl;
				std::vector<totempole_sorted_entry_type> entries;
				if(!this->toi_reader.findOverlap(interval.contigID, interval.start, interval.stop, entries)){
					//std::cerr << "could not find: " << interval << std::endl;
					continue;
				}

				for(U32 i = 0; i < entries.size(); ++i){
					//std::cerr << "contigID:pos-pos found: " << entries[i] << std::endl;
					this->interval_totempole_enties->push_back(entries[i]);
				}
			}
		}
	}

	return(this->interval_totempole_enties->size() > 0);
}

bool TomahawkOutputReader::__Open(const std::string input){
	this->stream_.open(input, std::ios::binary | std::ios::in | std::ios::ate);
	if(!this->stream_.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to open file: " << input << std::endl;
		return false;
	}

	this->filesize = this->stream_.tellg();
	this->stream_.seekg(0);

	if(!this->stream_.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Bad stream!" << std::endl;
		return false;
	}

	this->stream_ >> this->header;
	if(!this->header.validate(Tomahawk::Constants::WRITE_HEADER_LD_MAGIC)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to validate header!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::Open(const std::string input){
	if(!this->__Open(input))
		return false;

	if(!this->ParseHeader()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to parse header!" << std::endl;
		return false;
	}

	if(this->toi_reader.Open(input + "." + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX, this->contigs)){
		this->hasIndex = true;
	}

	return true;
}

bool TomahawkOutputReader::OpenExtend(const std::string input){
	if(!this->__Open(input))
		return false;

	if(!this->ParseHeaderExtend()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to extend header!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::__concat(const std::vector<std::string>& files, const std::string& output){
	if(files.size() == 0){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "No input files..." << std::endl;
		return false;
	}

	// open first one
	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "CONCAT") << "Opening input: " << files[0] << "..." << std::endl;

	if(!this->Open(files[0])){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Failed to parse: " << files[0] << "..." << std::endl;
		return false;
	}

	this->setWriterType(0);
	this->setWriteHeader(true);
	this->literals += "\n##tomahawk_concatCommand=" + Helpers::program_string();
	this->literals += "\n##tomahawk_concatFiles=";
	for(U32 i = 0; i < files.size(); ++i)
		this->literals += files[i] + ',';

	if(!this->OpenWriter(output)){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Failed to open writer..." << std::endl;
		return false;
	}

	while(this->parseBlock()){
		this->writer->write(this->data_);
	}

	for(U32 i = 1; i < files.size(); ++i){
		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "CONCAT") << "Opening input: " << files[i] << "..." << std::endl;

		this->stream_.close();
		if(!this->OpenExtend(files[i])){
			std::cerr << Helpers::timestamp("ERROR","TWO") << "Failed to parse: " << files[i] << "..." << std::endl;
			return false;
		}

		while(this->parseBlock()){
			this->writer->write(this->data_);
		}
	}

	this->writer->flush();
	this->writer->close();
	return true;
}

bool TomahawkOutputReader::concat(const std::vector<std::string>& files, const std::string& output){
	if(files.size() == 0){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "No input files given..." << std::endl;
		return false;
	}

	return(this->__concat(files, output));
}

bool TomahawkOutputReader::concat(const std::string& file_list, const std::string& output){
	if(file_list.size() == 0){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "No input file list given..." << std::endl;
		return false;
	}

	std::ifstream file_list_read(file_list);
	if(!file_list_read.good()){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Failed to get file_list..." << std::endl;
		return false;
	}

	std::vector<std::string> files;
	std::string line;
	while(getline(file_list_read, line)){
		if(line.size() == 0){
			std::cerr << Helpers::timestamp("WARNING","TWO") << "Empty line" << std::endl;
			break;
		}
		files.push_back(line);
	}

	return(this->__concat(files, output));
}

bool TomahawkOutputReader::ParseHeader(void){
	if(this->header.n_contig == 0)
		return false;

	if(this->header.n_contig < 1024)
		this->contig_htable = new hash_table(1024);
	else
		this->contig_htable = new hash_table(this->header.n_contig * 2);

	this->contigs = new contig_type[this->header.n_contig];
	U32* ret;

	for(U32 i = 0; i < this->header.n_contig; ++i){
		this->stream_ >> this->contigs[i];
		if(!this->contig_htable->GetItem(&this->contigs[i].name[0], &this->contigs[i].name, ret, this->contigs[i].name.size())){
			// Add to hash table
			this->contig_htable->SetItem(&this->contigs[i].name[0], &this->contigs[i].name, i, this->contigs[i].name.size());
		} else {
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Duplicated contig name: " << this->contigs[i].name << "!" << std::endl;
			exit(1); // unrecoverable error
		}
	}

	if(!this->tgzf_controller.InflateBlock(this->stream_, this->compressed_buffer)){
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "Failed to get TWO block" << std::endl;
		return false;
	}

	this->literals = std::string(this->tgzf_controller.buffer.data());

	return true;
}

bool TomahawkOutputReader::ParseHeaderExtend(void){
	if(this->header.n_contig == 0)
		return false;

	U32* ret;
	for(U32 i = 0; i < this->header.n_contig; ++i){
		this->stream_ >> this->contigs[i];
		// std::cerr << this->contigs[i] << std::endl;
		if(!this->contig_htable->GetItem(&this->contigs[i].name[0], &this->contigs[i].name, ret, this->contigs[i].name.size())){
			std::cerr << Helpers::timestamp("ERROR","TWO") << "Contig does not exist in other file" << std::endl;
			return false;
		}
	}

	if(!this->tgzf_controller.InflateBlock(this->stream_, this->compressed_buffer)){
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "Failed to get TWO block" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::seekBlock(const U32 blockID){
	if(this->toi_reader.ERROR_STATE != toi_reader_type::TOI_OK){
		std::cerr << Helpers::timestamp("ERROR","TOI") << "Index is bad!" << std::endl;
		return false;
	}

	if(blockID > this->toi_reader.size()){
		std::cerr << Helpers::timestamp("ERROR","TOI") << "Illegal blockID (" << blockID << ">" << this->toi_reader.size() << ")!" << std::endl;
		return false;
	}

	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Stream is bad!" << std::endl;
		return false;
	}

	this->stream_.seekg(this->toi_reader[blockID].getStartOffset());
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Stream is bad following seek!" << std::endl;
		return false;
	}

	return(this->parseBlock());
}

bool TomahawkOutputReader::summary(const std::string& input, const U32 bins){
	TWO::TomahawkOutputStatsContainer container(bins);

	// Natural output required parsing
	while(this->parseBlock()){
		output_container_reference_type o(this->data_);
		for(U32 i = 0; i < o.size(); ++i)
			container += o[i];
	}

	std::cerr << "R2\t" << container.R2.within.getTotal() << '\t' << container.R2.across.getTotal() << '\t' << container.R2.global.getTotal()  << std::endl;
	std::cerr << container.R2 << std::endl;
	std::cerr << "D\t" << container.D.within.getTotal() << '\t' << container.D.across.getTotal() << '\t' << container.D.global.getTotal()  << std::endl;
	std::cerr << container.D << std::endl;
	std::cerr << "Dprime\t" << container.Dprime.within.getTotal() << '\t' << container.Dprime.across.getTotal() << '\t' << container.Dprime.global.getTotal()  << std::endl;
	std::cerr << container.Dprime << std::endl;

	return true;
}

bool TomahawkOutputReader::index(const std::string& input){
	std::vector<std::string> paths = Helpers::filePathBaseExtension(input);
	std::string basePath = paths[0];
	std::string baseName;

	if(basePath.size() > 0)
		basePath += '/';

	if(paths[3].size() == Tomahawk::Constants::OUTPUT_LD_SUFFIX.size() &&
	   strncasecmp(&paths[3][0], &Tomahawk::Constants::OUTPUT_LD_SUFFIX[0], Tomahawk::Constants::OUTPUT_LD_SUFFIX.size()) == 0)
		baseName = paths[2];
	else baseName = paths[1];

	// Open writer
	// Set controller
	toi_header_type toi_header(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, this->header.samples, this->header.n_contig);
	// We assume data is expanded and sorted
	toi_header.controller.sorted = 1;
	toi_header.controller.expanded = 1;
	toi_header.controller.partial_sort = 0;

	twoi_writer_type writer(this->contigs, &this->header, toi_header);
	writer.open(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX);


	return true;
}

bool TomahawkOutputReader::setWriterType(const int type){
	if(type == 0)
		this->writer_output_type = WRITER_TYPE::binary;
	else if(type == 1)
		this->writer_output_type = WRITER_TYPE::natural;
	else {
		std::cerr << Tomahawk::Helpers::timestamp("ERROR","READER") << "Unknown writer type: " << type << std::endl;
		return false;
	}
	return true;
}
*/

} /* namespace Tomahawk */
