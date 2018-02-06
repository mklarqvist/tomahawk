#ifndef TOMAHAWK_TOMAHAWKREADER_H_
#define TOMAHAWK_TOMAHAWKREADER_H_

#include <string>
#include <fstream>
#include <algorithm>
#include <bitset>

#include "../support/MagicConstants.h"
#include "../io/compression/TGZFController.h"
#include "../io/compression/GZFConstants.h"
#include "../interface/Timer.h"
#include "../interface/ProgressBar.h"
#include "../algorithm/Balancer.h"
#include "base/meta_entry.h"
#include "base/twk_reader_implementation.h"
#include "ld_calculation_slave.h"
#include "TomahawkCalcParameters.h"

namespace Tomahawk {

// TomahawkReader class simply reads compressed data from disk
class TomahawkReader {
	typedef TomahawkCalcParameters    parameter_type;
	typedef Totempole::TotempoleEntry totempole_entry;

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

	bool Open(const std::string input);

	// Reader functions
	bool getBlocks(void);
	bool getBlocks(std::vector<U32>& blocks);
	bool getBlocks(std::vector< std::pair<U32, U32> >& blocks);
	bool getBlock(const U32 blockID);

	// Output functions
	bool outputBlocks(std::vector<U32>& blocks);
	bool outputBlocks();

	inline const BYTE& getBitWidth(void) const{ return(this->bit_width_); }
	inline Totempole::TotempoleReader& getTotempole(void){ return(this->totempole_); }
	inline const DataOffsetPair& getOffsetPair(const U32 p) const{ return(this->blockDataOffsets_[p]); }
	inline const size_t DataOffsetSize(void) const{ return(this->blockDataOffsets_.size()); }
	inline void setDropGenotypes(const bool yes){ this->dropGenotypes = yes; }
	inline void setShowHeader(const bool yes){ this->showHeader = yes; }

private:
	void DetermineBitWidth(void);
	template <class T> bool outputBlock(const U32 blockID);
	template <class T> bool WriteBlock(const char* data, const U32 blockID);
	bool Validate(void);
	bool ValidateHeader(std::ifstream& in) const;

private:
	U64 samples;     // has to match header
	float version;   // has to match header
	U64 filesize_;   // filesize
	BYTE bit_width_; // bit width
	bool dropGenotypes;
	bool showHeader; // flag to output header or not
	std::ifstream stream_; // reader stream
	Totempole::TotempoleReader totempole_;


	IO::BasicBuffer buffer_; // input buffer
	IO::BasicBuffer data_; // inflate buffer
	IO::BasicBuffer outputBuffer_; // output buffer
	IO::TGZFController tgzf_controller_;

	std::vector<DataOffsetPair> blockDataOffsets_;

	IO::GenericWriterInterace* writer;
};

template <class T>
bool TomahawkReader::outputBlock(const U32 blockID){
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Stream bad " << blockID << std::endl;
		return false;
	}

	//std::cerr << "getblock " << blockID  << " seek to " << this->totempole_[blockID].byte_offset << std::endl;
	this->stream_.seekg(this->totempole_[blockID].byte_offset);
	if(!this->stream_.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Failed search..." << std::endl;
		return false;
	}

	// Determine byte-width of data
	U32 readLength = 0;
	if(blockID != this->totempole_.getBlocks() - 1)
		readLength = this->totempole_[blockID + 1].byte_offset - this->totempole_[blockID].byte_offset;
	else
		readLength = this->filesize_ - Constants::eof_length*sizeof(U64) - this->totempole_[this->totempole_.getBlocks()-1].byte_offset;

	if(readLength > this->buffer_.capacity()){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Impossible: " << readLength << '/' << this->buffer_.capacity() << std::endl;
		exit(1);
	}

	// Read from start to start + byte-width
	if(!this->stream_.read(&this->buffer_.data[0], readLength)){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Failed read: " << this->stream_.good() << '\t' << this->stream_.fail() << '/' << this->stream_.eof() << std::endl;
		std::cerr << this->stream_.gcount() << '/' << readLength << std::endl;
		return false;
	}
	// Set buffer byte-width to data loaded
	this->buffer_.pointer = readLength;

	// Keep track of position because inflate function moves pointer
	char* data_position = &this->data_.data[this->data_.pointer];

	// Inflate TGZF block
	if(!this->tgzf_controller_.Inflate(this->buffer_, this->data_)){
		std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to inflate DATA..." << std::endl;
		return false;
	}

	// Todo: move to function
	this->WriteBlock<T>(data_position, blockID);

	return true;
}

template <class T>
bool TomahawkReader::WriteBlock(const char* data, const U32 blockID){
	//Todo
	/*
	TomahawkBlock<T> tomahawk_controller(data, this->totempole_[blockID]);

	// For each variant in Tomahawk block
	for(U32 j = 0; j < tomahawk_controller.support->variants; ++j){
		tomahawk_controller.WriteVariant(this->totempole_, this->outputBuffer_, this->dropGenotypes);

		// Next variant
		++tomahawk_controller;

		// Keep flushing regularly

		if(this->outputBuffer_.size() > 65536){
			//this->writer->write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
			std::cout.write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
			this->outputBuffer_.reset();
		}
	}

	// Flush last
	//this->writer->write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
	std::cout.write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);

	// Reset buffers
	this->outputBuffer_.reset(); // reset
	this->data_.reset(); // reset
	*/

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKREADER_H_ */
