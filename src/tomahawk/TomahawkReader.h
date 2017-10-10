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
#include "TomahawkCalcParameters.h"
#include "base/TomahawkEntryMeta.h"
#include "TomahawkCalculateSlave.h"

namespace Tomahawk {

// TomahawkReader class simply reads compressed data from disk
class TomahawkReader {
	typedef TomahawkCalcParameters parameter_type;
	typedef Totempole::TotempoleEntry totempole_entry;
	typedef IO::TGZFController tgzf_controller_type;
	typedef Totempole::TotempoleReader totempole_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::GenericWriterInterace writer_interface;

public:
	// Used to keep track of char pointer offsets in buffer
	// and what Totempole entry is associated with that position
	// given the data that is loaded
	// This is equivalent to Totempole virtual offsets given
	// the data loaded in memory
	struct DataOffsetPair{
		DataOffsetPair(const char* data, const totempole_entry& entry) : entry(entry), data(data){}
		~DataOffsetPair(){}

		const totempole_entry& entry;
		const char* data;
	};

private:
	typedef std::vector<DataOffsetPair> offset_vector;

public:
	TomahawkReader();
	virtual ~TomahawkReader();

	bool Open(const std::string input);

	// Reader functions
	bool nextBlock(const bool clear = true);
	bool getBlocks(void);
	bool getBlocks(std::vector<U32>& blocks);
	bool getBlocks(std::vector< std::pair<U32, U32> >& blocks);
	bool getBlock(const U32 blockID);

	// Output functions
	bool outputBlocks(std::vector<U32>& blocks);
	bool outputBlocks();

	inline const BYTE& getBitWidth(void) const{ return(this->bit_width); }
	inline Totempole::TotempoleReader& getTotempole(void){ return(this->totempole); }
	inline const DataOffsetPair& getOffsetPair(const U32 p) const{ return(this->blockDataOffsets[p]); }
	inline const size_t DataOffsetSize(void) const{ return(this->blockDataOffsets.size()); }
	inline void setDropGenotypes(const bool yes){ this->dropGenotypes = yes; }
	inline void setShowHeader(const bool yes){ this->showHeader = yes; }

private:
	void DetermineBitWidth(void);
	template <class T> bool outputBlock(const U32 blockID);
	template <class T> bool WriteBlock(const char* data, const U32 blockID);
	bool Validate(void);
	bool ValidateHeader(std::ifstream& in) const;
	template <class T> bool __calculateTajimaD(const U32 bin_size);
	template <class T> bool __calculateFST(void);
	template <class T> bool __calculateSFS(void);

protected:
	U64 samples;     // has to match header
	float version;   // has to match header
	U64 filesize;   // filesize
	BYTE bit_width; // bit width
	bool dropGenotypes; // drop genotypes in view mode
	bool showHeader; // flag to output header or not

	U32 currentBlockID; // for iterator

	std::ifstream stream; // reader stream

	totempole_type totempole; // totempole reader
	buffer_type buffer; // input buffer
	buffer_type data; // inflate buffer
	buffer_type outputBuffer; // output buffer
	tgzf_controller_type tgzf_controller; // tgzf controller
	offset_vector blockDataOffsets; // internal virtual offsets into buffer
	writer_interface* writer; // writer interface
};

template <class T>
bool TomahawkReader::outputBlock(const U32 blockID){
	if(!this->stream.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Stream bad " << blockID << std::endl;
		return false;
	}

	//std::cerr << "getblock " << blockID  << " seek to " << this->totempole_[blockID].byte_offset << std::endl;
	this->stream.seekg(this->totempole[blockID].byte_offset);
	if(!this->stream.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Failed search..." << std::endl;
		return false;
	}

	// Determine byte-width of data
	const U32 readLength = this->totempole[blockID].byte_offset_end - this->totempole[blockID].byte_offset;

	if(readLength > this->buffer.capacity()){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Impossible: " << readLength << '/' << this->buffer.capacity() << std::endl;
		exit(1);
	}

	// Read from start to start + byte-width
	if(!this->stream.read(&this->buffer.data[0], readLength)){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Failed read: " << this->stream.good() << '\t' << this->stream.fail() << '/' << this->stream.eof() << std::endl;
		std::cerr << this->stream.gcount() << '/' << readLength << std::endl;
		return false;
	}
	// Set buffer byte-width to data loaded
	this->buffer.pointer = readLength;

	// Keep track of position because inflate function moves pointer
	char* data_position = &this->data.data[this->data.pointer];

	// Inflate TGZF block
	if(!this->tgzf_controller.Inflate(this->buffer, this->data)){
		std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to inflate DATA..." << std::endl;
		return false;
	}

	// Todo: move to function
	this->WriteBlock<T>(data_position, blockID);

	return true;
}

template <class T>
bool TomahawkReader::WriteBlock(const char* data, const U32 blockID){
	TomahawkIterator<T> tomahawk_controller(data, this->totempole[blockID]);

	// For each variant in Tomahawk block
	for(U32 j = 0; j < tomahawk_controller.support->variants; ++j){
		tomahawk_controller.WriteVariant(this->totempole, this->outputBuffer, this->dropGenotypes);

		// Next variant
		++tomahawk_controller;

		// Keep flushing regularly
		// arbitrary threshold at 65536 bytes
		if(this->outputBuffer.size() > 65536){
			//this->writer->write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
			std::cout.write(&this->outputBuffer.data[0], this->outputBuffer.pointer);
			this->outputBuffer.reset();
		}
	}

	// Flush last
	//this->writer->write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
	std::cout.write(&this->outputBuffer.data[0], this->outputBuffer.pointer);

	// Reset buffers
	this->outputBuffer.reset(); // reset
	this->data.reset(); // reset

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKREADER_H_ */
