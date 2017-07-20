#ifndef TOMAHAWK_TOMAHAWKREADER_H_
#define TOMAHAWK_TOMAHAWKREADER_H_

#include <string>
#include <fstream>
#include <algorithm>
#include <bitset>

#include "MagicConstants.h"
#include "../io/GZController.h"
#include "TomahawkEntryMeta.h"
#include "../totempole/TotempoleReader.h"
#include "../io/IOConstants.h"
#include "TomahawkCalculateSlave.h"
#include "../interface/Timer.h"
#include "../interface/ProgressBar.h"
#include "TomahawkCalcParameters.h"
#include "../algorithm/Balancer.h"

// MAJOR TODO: SEPARATE OUT TWK READER FROM CALC FUNCTIONS

namespace Tomahawk {

// TomahawkReader class simply reads compressed data from disk
class TomahawkReader {
	typedef TomahawkCalcParameters parameter_type;
	typedef TotempoleEntry totempole_entry;

public:
	// Used to keep track of char pointer offsets in buffer
	// and what Totempole entry is associated with that position
	struct DataOffsetPair{
		DataOffsetPair(const char* data, const totempole_entry& entry) : entry(entry), data(data){}
		~DataOffsetPair(){}

		const totempole_entry& entry;
		const char* data;
	};

public:
	TomahawkReader();
	~TomahawkReader();

	// Reader functions
	bool Open(const std::string input);
	bool getBlocks(void);
	bool getBlocks(std::vector<U32>& blocks);
	bool getBlocks(std::vector< std::pair<U32, U32> >& blocks);
	bool getBlock(const U32 blockID);
	bool outputBlocks(std::vector<U32>& blocks);
	bool outputBlocks();

	// Todo: reimplement for writing Tomahawk output data
	//bool SelectWriterOutputType(const IO::GenericWriterInterace::type writer_type);
	//void SetOutputType(IO::GenericWriterInterace::compression type){ this->parameters.compression_type = type; }
	//bool OpenWriter(void);
	//bool OpenWriter(const std::string destination);

	inline const BYTE& getBitWidth(void) const{ return(this->bit_width_); }
	inline const TotempoleReader& getTotempole(void) const{ return(this->totempole_); }
	inline const DataOffsetPair& getOffsetPair(const U32 p) const{ return(this->blockDataOffsets_[p]); }
	inline const size_t DataOffsetSize(void) const{ return(this->blockDataOffsets_.size()); }

private:
	void DetermineBitWidth(void);
	U64 GetUncompressedSizes(std::vector<U32>& blocks);
	U64 GetUncompressedSizes(std::vector< std::pair<U32, U32> >& blocks);
	U64 GetUncompressedSizes(void);
	template <class T> bool outputBlock(const U32 blockID);

	// Reader functions
	template <class T> bool WriteBlock(const char* data, const U32 blockID);
	bool Validate(void);
	bool ValidateHeader(std::ifstream& in) const;

private:
	U64 samples; // has to match header
	float version; // has to match header
	U64 filesize_;
	BYTE bit_width_;
	std::ifstream stream_; // reader stream

	IO::GenericWriterInterace* writer;
	TotempoleReader totempole_;

	// Todo: Move out
	IO::BasicBuffer buffer_;
	IO::BasicBuffer data_;
	IO::BasicBuffer outputBuffer_;
	IO::GZController gzip_controller_;

	std::vector<DataOffsetPair> blockDataOffsets_;
};

template <class T>
bool TomahawkReader::outputBlock(const U32 blockID){
	if(!this->stream_.good()){
		std::cerr << "stream bad " << blockID << std::endl;
		return false;
	}

	//std::cerr << "getblock " << blockID  << " seek to " << this->totempole_[blockID].byte_offset << std::endl;
	this->stream_.seekg(this->totempole_[blockID].byte_offset);
	if(!this->stream_.good()){
		std::cerr << "Failed search" << std::endl;
		return false;
	}

	// Determine byte-width of data
	U32 readLength = 0;
	if(blockID != this->totempole_.getBlocks() - 1)
		readLength = this->totempole_[blockID + 1].byte_offset - this->totempole_[blockID].byte_offset;
	else
		readLength = this->filesize_ - Constants::eof_length*sizeof(U64) - this->totempole_[this->totempole_.getBlocks()-1].byte_offset;

	if(readLength > this->buffer_.capacity()){
		std::cerr << "impossible: " << readLength << '/' << this->buffer_.capacity() << std::endl;
		exit(1);
	}

	// Read from start to start + byte-width
	if(!this->stream_.read(&this->buffer_.data[0], readLength)){
		std::cerr << "Failed read: " << this->stream_.good() << '\t' << this->stream_.fail() << '/' << this->stream_.eof() << std::endl;
		std::cerr << this->stream_.gcount() << '/' << readLength << std::endl;
		return false;
	}
	// Set buffer byte-width to data loaded
	this->buffer_.pointer = readLength;

	// Keep track of position because inflate function moves pointer
	char* data_position = &this->data_.data[this->data_.pointer];

	// Inflata TGZF block
	if(!this->gzip_controller_.Inflate(this->buffer_, this->data_)){
		std::cerr << "failed" << std::endl;
		return false;
	}

	// Todo: move to function
	this->WriteBlock<T>(data_position, blockID);

	return true;
}

template <class T>
bool TomahawkReader::WriteBlock(const char* data, const U32 blockID){
	TomahawkBlock<T> tomahawk_controller(data, this->totempole_[blockID]);

	// For each variant in Tomahawk block
	for(U32 j = 0; j < tomahawk_controller.support->variants; ++j){
		tomahawk_controller.WriteVariant(this->totempole_, this->outputBuffer_);

		// Next variant
		++tomahawk_controller;

		// Keep flushing regularly
		if(this->outputBuffer_.size() > 65536){
			this->writer->write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
			//std::cout.write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
			this->outputBuffer_.reset();
		}
	}

	// Flush last
	this->writer->write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);

	// Reset buffers
	this->outputBuffer_.reset(); // reset
	this->data_.reset(); // reset

	return true;
}



} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKREADER_H_ */
