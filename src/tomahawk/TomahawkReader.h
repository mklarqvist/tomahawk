#ifndef TOMAHAWK_TOMAHAWKREADER_H_
#define TOMAHAWK_TOMAHAWKREADER_H_

#include <string>
#include <fstream>
#include <algorithm>
#include <bitset>
#include <regex>

#include "../interface/progressbar.h"
#include "../support/MagicConstants.h"
#include "../io/compression/TGZFController.h"
#include "../io/compression/GZFConstants.h"
#include "../interface/timer.h"
#include "meta_entry.h"
#include "twk_reader_implementation.h"
#include "ld_calculation_slave.h"
#include "TomahawkCalcParameters.h"
#include "../index/index.h"
#include "genotype_container_reference.h"
#include "../index/tomahawk_header.h"
#include "../third_party/intervalTree.h"

namespace Tomahawk {

// TomahawkReader class simply reads compressed data from disk
class TomahawkReader {
	typedef TomahawkCalcParameters     parameter_type;
	typedef Totempole::IndexEntry      totempole_entry;
	typedef TomahawkHeader             header_type;
	typedef Index                      index_type;
	typedef IO::BasicBuffer            buffer_type;
	typedef IO::TGZFController         tgzf_controller_type;
	typedef Totempole::Footer          footer_type;
	typedef Algorithm::ContigInterval  interval_type;
	typedef Algorithm::IntervalTree<interval_type, U32> tree_type;


public:
	// Used to keep track of char pointer offsets in buffer
	// and what Totempole entry is associated with that position
	struct DataOffsetPair{
		DataOffsetPair(const char* data, const U64 l_buffer, const totempole_entry& entry) : entry(entry), l_buffer(l_buffer), data(data){}
		~DataOffsetPair(){}

		const totempole_entry& entry;
		const U64 l_buffer;
		const char* data;
	};

public:
	TomahawkReader();
	~TomahawkReader();

	bool open(const std::string input);

	// Reader functions
	bool getBlocks(void);
	bool getBlocks(std::vector<U32>& blocks);
	bool getBlocks(std::vector< std::pair<U32, U32> >& blocks);
	bool getBlock(const U32 blockID);

	// Accessors
	inline footer_type& getFooter(void){ return(this->footer_); }
    inline const footer_type& getFooter(void) const{ return(this->footer_); }
	inline const index_type& getIndex(void) const{ return(*this->index_); }
	inline index_type& getIndex(void){ return(*this->index_); }
	inline const header_type& getHeader(void) const{ return(this->header_); }
	inline header_type& getHeader(void){ return(this->header_); }
	inline index_type* getIndexPointer(void){ return(this->index_); }

	// Output functions
	bool outputBlocks(std::vector<U32>& blocks);
	bool outputBlocks();

	bool addRegions(const std::vector<std::string>& intervals){
		if(intervals.size() == 0)
			return true;

		if(this->interval_tree_entries == nullptr)
			this->interval_tree_entries = new std::vector<interval_type>[this->getHeader().getMagic().getNumberContigs()];

		if(this->interval_tree == nullptr){
			this->interval_tree = new tree_type*[this->getHeader().getMagic().getNumberContigs()];
			for(U32 i = 0; i < this->getHeader().getMagic().getNumberContigs(); ++i)
				this->interval_tree[i] = nullptr;
		}

		// Parse intervals
		for(U32 i = 0; i < intervals.size(); ++i){
			interval_type interval;
			// has colon (:)
			if(intervals[i].find(':') != std::string::npos){
				std::vector<std::string> first = Helpers::split(intervals[i], ':');
				if(first.size() != 2){
					std::cerr << Helpers::timestamp("ERROR") << "Malformed interval string: " << intervals[i] << "..." << std::endl;
					return false;
				}

				const S32 contigID = this->header_.getContigID(first[0]);
				if(contigID == -1){
					std::cerr << Helpers::timestamp("ERROR") << "Contig: " << intervals[i] << " is not defined in this file..." << std::endl;
					return false;
				}

				interval.state = Algorithm::ContigInterval::INTERVAL_FULL;
				interval.contigID = contigID;
				std::vector<std::string> sections = Helpers::split(first[1],'-');
				if(sections.size() == 2){ // Is contig + position->position
					if(std::regex_match(sections[0], std::regex("^[0-9]{1,}([\\.]{1}[0-9]{1,})?([eE]{1}[0-9]{1,})?$"))){
						interval.start = atof(sections[0].data()) - 1;
					} else {
						std::cerr << Helpers::timestamp("ERROR") << "Illegal number: " << sections[0] << " in interval string " << intervals[i] << std::endl;
						return false;
					}

					if(std::regex_match(sections[1], std::regex("^[0-9]{1,}([\\.]{1}[0-9]{1,})?([eE]{1}[0-9]{1,})?$"))){
						interval.stop = atof(sections[1].data());
					} else {
						std::cerr << Helpers::timestamp("ERROR") << "Illegal number: " << sections[1] << " in interval string " << intervals[i] << std::endl;
						return false;
					}

				} else if(sections.size() == 1){ // Is contig + position
					interval.state = Algorithm::ContigInterval::INTERVAL_POSITION;
					if(std::regex_match(sections[0], std::regex("^[0-9]{1,}([\\.]{1}[0-9]{1,})?([eE]{1}[0-9]{1,})?$"))){
						interval.start = atof(sections[0].data()) - 1;
						interval.stop  = atof(sections[0].data());
					} else {
						std::cerr << Helpers::timestamp("ERROR") << "Illegal number: " << sections[0] << " in interval string " << intervals[i] << std::endl;
						return false;
					}

				} else {
					std::cerr << Helpers::timestamp("ERROR") << "Malformed interval string: " << intervals[i] << std::endl;
					return false;
				}
			} else { // Is contig  only
				interval.state = Algorithm::ContigInterval::INTERVAL_CONTIG_ONLY;
				const S32 contigID = this->header_.getContigID(intervals[i]);
				if(contigID == -1){
					std::cerr << Helpers::timestamp("ERROR") << "Contig: " << intervals[i] << " is not defined in this file..." << std::endl;
					return false;
				}
				interval.contigID = contigID; // Todo
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


	inline const BYTE& getBitWidth(void) const{ return(this->bit_width_); }
	inline const DataOffsetPair& getOffsetPair(const U32 p) const{ return(this->blockDataOffsets_[p]); }
	inline const size_t DataOffsetSize(void) const{ return(this->blockDataOffsets_.size()); }

	inline void setDropGenotypes(const bool yes){ this->dropGenotypes = yes; }
	inline void setShowHeader(const bool yes){ this->showHeader = yes; }

private:
	void DetermineBitWidth(void);
	template <class T> bool outputBlock(const U32 blockID);
	template <class T> bool outputBlockFilter(const U32 blockID);
	template <class T> bool WriteBlock(const char* data, const U32 blockID);
	template <class T> bool WriteBlockFilter(const char* data, const U32 blockID);

private:
	U64            filesize_;  // filesize
	U64            offset_end_of_data_;
	BYTE           bit_width_; // bit width
	bool           dropGenotypes;
	bool           showHeader; // flag to output header or not
	std::ifstream  stream_;    // reader stream

	header_type header_;
	footer_type footer_;
	index_type* index_;

	buffer_type          buffer_;       // input buffer
	buffer_type          data_;         // inflate buffer
	buffer_type          outputBuffer_; // output buffer
	tgzf_controller_type tgzf_controller_;

	std::vector<DataOffsetPair> blockDataOffsets_;

	IO::GenericWriterInterace* writer;

	tree_type** interval_tree;
	std::vector<interval_type>* interval_tree_entries;
};

template <class T>
bool TomahawkReader::outputBlock(const U32 blockID){
	if(!this->getBlock(blockID)){
		std::cerr << "failed to get block" << std::endl;
		return false;
	}
	this->WriteBlock<T>(this->data_.data(), blockID);

	return true;
}

template <class T>
bool TomahawkReader::outputBlockFilter(const U32 blockID){
	if(!this->getBlock(blockID)){
		std::cerr << "failed to get block" << std::endl;
		return false;
	}
	this->WriteBlockFilter<T>(this->data_.data(), blockID);

	return true;
}

template <class T>
bool TomahawkReader::WriteBlockFilter(const char* const data, const U32 blockID){
	Base::GenotypeContainerReference<T> o(data,
                                          this->index_->getContainer()[blockID].uncompressed_size,
                                          this->index_->getContainer()[blockID],
                                          this->header_.magic_.getNumberSamples(),
                                          false);

	if(this->interval_tree == nullptr)
		return false;

	// For each variant in Tomahawk block
	for(U32 j = 0; j < o.size(); ++j){
		if(this->interval_tree[this->index_->getContainer()[blockID].contigID] == nullptr){
			++o;
			continue;
		}

		std::vector<interval_type> ret = this->interval_tree[this->index_->getContainer()[blockID].contigID]->findOverlapping(o.currentMeta().position, o.currentMeta().position);
		if(ret.size() == 0){
			++o;
			continue;
		}

		const char separator = o.currentMeta().getPhaseVCFCharacter();

		this->outputBuffer_ += this->header_.contigs_[this->index_->getContainer()[blockID].contigID].name;
		this->outputBuffer_ += '\t';
		this->outputBuffer_ += std::to_string(o.currentMeta().position);
		this->outputBuffer_ += "\t.\t";
		this->outputBuffer_ += o.currentMeta().getRefAllele();
		this->outputBuffer_ += '\t';
		this->outputBuffer_ += o.currentMeta().getAltAllele();
		this->outputBuffer_ += "\t.\t.\t";
		this->outputBuffer_ += "HWE_P=";
		this->outputBuffer_ += std::to_string(o.currentMeta().HWE_P);
		this->outputBuffer_ += ";MAF=";
		this->outputBuffer_ += std::to_string(o.currentMeta().AF);
		if(this->dropGenotypes == false){
			this->outputBuffer_ += "\tGT\t";
			for(U32 i = 0; i < o.currentMeta().runs; ++i){
				if(i != 0) this->outputBuffer_ += '\t';
				const char& left  = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[o[i].alleleA];
				const char& right = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[o[i].alleleB];
				this->outputBuffer_ += left;
				this->outputBuffer_ += separator;
				this->outputBuffer_ += right;
				for(U32 r = 1; r < o[i].runs; ++r){
					this->outputBuffer_ += '\t';
					this->outputBuffer_ += left;
					this->outputBuffer_ += separator;
					this->outputBuffer_ += right;
				}
			}
		}
		this->outputBuffer_ += '\n';
		++o;

		if(this->outputBuffer_.size() > 65536){
			//this->writer->write(this->outputBuffer_.data(), this->outputBuffer_.n_chars);
			std::cout.write(this->outputBuffer_.data(), this->outputBuffer_.n_chars);
			this->outputBuffer_.reset();
		}
	}

	// Flush last
	//this->writer->write(this->outputBuffer_.data(), this->outputBuffer_.n_chars);
	std::cout.write(this->outputBuffer_.data(), this->outputBuffer_.n_chars);

	// Reset buffers
	this->outputBuffer_.reset(); // reset
	this->data_.reset(); // reset

	return true;
}

template <class T>
bool TomahawkReader::WriteBlock(const char* const data, const U32 blockID){
	Base::GenotypeContainerReference<T> o(data,
                                          this->index_->getContainer()[blockID].uncompressed_size,
                                          this->index_->getContainer()[blockID],
                                          this->header_.magic_.getNumberSamples(),
                                          false);

	// For each variant in Tomahawk block
	for(U32 j = 0; j < o.size(); ++j){
		const char separator = o.currentMeta().getPhaseVCFCharacter();

		this->outputBuffer_ += this->header_.contigs_[this->index_->getContainer()[blockID].contigID].name;
		this->outputBuffer_ += '\t';
		this->outputBuffer_ += std::to_string(o.currentMeta().position);
		this->outputBuffer_ += "\t.\t";
		this->outputBuffer_ += o.currentMeta().getRefAllele();
		this->outputBuffer_ += '\t';
		this->outputBuffer_ += o.currentMeta().getAltAllele();
		this->outputBuffer_ += "\t.\t.\t";
		this->outputBuffer_ += "HWE_P=";
		this->outputBuffer_ += std::to_string(o.currentMeta().HWE_P);
		this->outputBuffer_ += ";MAF=";
		this->outputBuffer_ += std::to_string(o.currentMeta().AF);
		if(this->dropGenotypes == false){
			this->outputBuffer_ += "\tGT\t";
			for(U32 i = 0; i < o.currentMeta().runs; ++i){
				if(i != 0) this->outputBuffer_ += '\t';
				const char& left  = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[o[i].alleleA];
				const char& right = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[o[i].alleleB];
				this->outputBuffer_ += left;
				this->outputBuffer_ += separator;
				this->outputBuffer_ += right;
				for(U32 r = 1; r < o[i].runs; ++r){
					this->outputBuffer_ += '\t';
					this->outputBuffer_ += left;
					this->outputBuffer_ += separator;
					this->outputBuffer_ += right;
				}
			}
		}
		this->outputBuffer_ += '\n';
		++o;

		if(this->outputBuffer_.size() > 65536){
			//this->writer->write(this->outputBuffer_.data(), this->outputBuffer_.n_chars);
			std::cout.write(this->outputBuffer_.data(), this->outputBuffer_.n_chars);
			this->outputBuffer_.reset();
		}
	}

	// Flush last
	//this->writer->write(this->outputBuffer_.data(), this->outputBuffer_.n_chars);
	std::cout.write(this->outputBuffer_.data(), this->outputBuffer_.n_chars);

	// Reset buffers
	this->outputBuffer_.reset(); // reset
	this->data_.reset(); // reset

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKREADER_H_ */
