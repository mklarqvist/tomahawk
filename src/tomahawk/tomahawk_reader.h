#ifndef TOMAHAWK_TOMAHAWK_READER_H_
#define TOMAHAWK_TOMAHAWK_READER_H_

#include <string>
#include <fstream>
#include <algorithm>
#include <bitset>
#include <regex>

#include "../support/magic_constants.h"
#include "io/compression/gz_constants.h"
#include "io/compression/tgzf_controller.h"
#include "interface/progressbar.h"
#include "interface/timer.h"
#include "meta_entry.h"
#include "twk_reader_implementation.h"
#include "ld_calculation_slave.h"
#include "index/index.h"
#include "genotype_container_reference.h"
#include "index/tomahawk_header.h"
#include "third_party/intervalTree.h"
#include "tomahawk_calc_parameters.h"

namespace tomahawk {

struct temp{
	temp() :
		n_hetA(0), n_hetB(0),
		n_homA(0), n_homB(0),
		n_missA(0), n_missB(0)
	{

	}

	void clear(void){
		this->n_hetA = 0; this->n_hetB = 0;
		this->n_homA = 0; this->n_homB = 0;
		this->n_missA = 0; this->n_missB = 0;
	}

	friend std::ostream& operator<<(std::ostream& os, const temp& self){
		os << self.n_hetA << "\t" << self.n_hetB << "\t" << self.n_homA << "\t" << self.n_homB << "\t" << self.n_missA << "\t" << self.n_missB;
		return(os);
	}

	U32 n_hetA, n_hetB;
	U32 n_homA, n_homB;
	U32 n_missA, n_missB;
};

// TomahawkReader class simply reads compressed data from disk
class TomahawkReader {
	typedef TomahawkCalcParameters     parameter_type;
	typedef totempole::IndexEntry      totempole_entry;
	typedef TomahawkHeader             header_type;
	typedef Index                      index_type;
	typedef io::BasicBuffer            buffer_type;
	typedef io::TGZFController         tgzf_controller_type;
	typedef totempole::Footer          footer_type;
	typedef algorithm::ContigInterval  interval_type;
	typedef algorithm::IntervalTree<interval_type, U32> tree_type;


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
	bool addRegions(const std::vector<std::string>& intervals);
	bool printHeader(std::ostream& stream) const;

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

	bool summaryIndividuals();

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

	template <class T> bool statsIndividual(std::vector<temp>& stats, const char* data, const U32 blockID);

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

	io::GenericWriterInterace* writer;

public:
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
	base::GenotypeContainerReference<T> o(data,
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
		this->outputBuffer_ += ";AF=";
		this->outputBuffer_ += std::to_string(o.currentMeta().AF);
		if(this->dropGenotypes == false){
			this->outputBuffer_ += "\tGT\t";
			for(U32 i = 0; i < o.currentMeta().runs; ++i){
				if(i != 0) this->outputBuffer_ += '\t';
				const char& left  = constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[o[i].alleleA];
				const char& right = constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[o[i].alleleB];
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
	base::GenotypeContainerReference<T> o(data,
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
		this->outputBuffer_ += ";AF=";
		this->outputBuffer_ += std::to_string(o.currentMeta().AF);
		if(this->dropGenotypes == false){
			this->outputBuffer_ += "\tGT\t";
			for(U32 i = 0; i < o.currentMeta().runs; ++i){
				if(i != 0) this->outputBuffer_ += '\t';
				const char& left  = constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[o[i].alleleA];
				const char& right = constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[o[i].alleleB];
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
bool TomahawkReader::statsIndividual(std::vector<temp>& stats, const char* const data, const U32 blockID){
	base::GenotypeContainerReference<T> o(data,
                                          this->index_->getContainer()[blockID].uncompressed_size,
                                          this->index_->getContainer()[blockID],
                                          this->header_.magic_.getNumberSamples(),
                                          false);



	//std::vector<temp> stats(this->getHeader().getMagic().n_samples);
	for(U32 j = 0; j < o.size(); ++j){
		U32 cumsum = 0;
		for(U32 i = 0; i < o.currentMeta().runs; ++i){
			for(U32 r = 0; r < o[i].runs; ++r, cumsum++){
				stats[cumsum].n_hetA += o[i].alleleA == 1;
				stats[cumsum].n_hetB += o[i].alleleB == 1;
				stats[cumsum].n_homA += o[i].alleleA == 0;
				stats[cumsum].n_homB += o[i].alleleB == 0;
				stats[cumsum].n_missA += o[i].alleleA == 2;
				stats[cumsum].n_missB += o[i].alleleB == 2;
			}
		}
		//std::cerr << cumsum << "/" << this->getHeader().getMagic().n_samples << std::endl;
		assert(cumsum == this->getHeader().getMagic().n_samples);
	}

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWK_READER_H_ */
