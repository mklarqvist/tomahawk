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

	// Stats
	bool calculateTajimaD(void);
	bool calculateFST(void);
	bool calculateSFS(void);
	bool calculateIBS(void);
	bool calculateROH(void);
	bool calculateNucleotideDiversity(void);

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
	template <class T> bool __calculateTajimaD(void);
	template <class T> bool __calculateFST(void);
	template <class T> bool __calculateSFS(void);

private:
	U64 samples;     // has to match header
	float version;   // has to match header
	U64 filesize;   // filesize
	BYTE bit_width; // bit width
	bool dropGenotypes;
	bool showHeader; // flag to output header or not

	U32 currentBlockID;

	std::ifstream stream; // reader stream
	Totempole::TotempoleReader totempole;


	IO::BasicBuffer buffer; // input buffer
	IO::BasicBuffer data; // inflate buffer
	IO::BasicBuffer outputBuffer; // output buffer
	IO::TGZFController tgzf_controller;

	std::vector<DataOffsetPair> blockDataOffsets;

	IO::GenericWriterInterace* writer;
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

template <class T>
bool TomahawkReader::__calculateTajimaD(void){
	this->buffer.resize(this->totempole.getLargestBlockSize() + 1);

	// Constants
	const double n = this->totempole.header.samples;
	double a1 = 0, a2 = 0;
	for(U32 i = 1; i < n - 1; ++i){
		a1 += (double)1/i;
		a2 += (double)1/(i*i);
	}

	const double b1 = (n + 1) / (3 * (n - 1));
	const double b2 = (2*(n*n + n + 3)) / (9*n * (n - 1));
	const double c1 = b1 - 1/a1;
	const double c2 = b2 - ((n + 2) / (a1*n)) + a2/(a1*a1);
	const double e1 = c1/a1;
	const double e2 = c2 / (a1*a1 + a2);

	// Variables
	T* lookup = new T[16];
	U32 counts = 0;
	U32 prevPos = 0;
	U64 f0 = 0;
	U64 f1 = 0;
	double comps = 0;

	U32 current_bin = 0;
	S32 previous_contigID = 0;
	U64 cumSum = 0;

	std::cout << "n_snps\tbinFrom\tbinTo\tD\tpi" << std::endl;
	for(U32 i = 0; i < this->totempole.header.blocks; ++i){
		if(!this->nextBlock()){
			std::cerr << "failed to get next block" << std::endl;
			return false;
		}

		// Initiate array
		memset(lookup, 0, sizeof(T)*16);

		// Now have data
		TomahawkIterator<T> controller(this->data.data, this->totempole[i]);
		const Support::TomahawkRunPacked<T>* runs = nullptr;
		const TomahawkEntryMeta<T>* meta = nullptr;

		while(controller.nextVariant(runs, meta)){
			for(U32 i = 0; i < meta->runs; ++i){
				lookup[runs[i].alleles] += runs[i].runs;
			}

			if((meta->position/1000)*1000 != current_bin || this->totempole[i].contigID != previous_contigID){
				if(counts > 0){
					const double S = counts;
					const double D = ((double)n/(n-1)*comps - (S/a1)) / sqrt(e1*S + e2*S*(S-1));
					std::cout << counts << '\t' << previous_contigID << '\t' << current_bin << '\t' << current_bin + 1000 << '\t' << cumSum << '\t' << cumSum+1000 << '\t' << D << '\t' << (double)n/(n-1)*comps << std::endl;
					cumSum += 1000;
				}

				// Reset
				counts = 0;
				f0 = 0;
				f1 = 0;
				comps = 0;
				current_bin = (meta->position/1000)*1000;
				previous_contigID = this->totempole[i].contigID;
			}

			f0 = 2*lookup[0] + lookup[1] + lookup[4];
			f1 = 2*lookup[5] + lookup[1] + lookup[4];
			comps += (double)f0*f1 / (((f0+f1) * (f0+f1) - (f0+f1)) / 2);
			++counts;
			prevPos = meta->position;

			memset(lookup, 0, sizeof(T)*16);
		}
	}

	// Last
	if(counts > 0){
		const double S = counts;
		const double D = ((double)n/(n-1)*comps - (S/a1)) / sqrt(e1*S + e2*S*(S-1));
		std::cout << counts << '\t' << previous_contigID << '\t' << current_bin << '\t' << current_bin + 1000 << '\t' << cumSum << '\t' << cumSum+1000 << '\t' << D << '\t' << (double)n/(n-1)*comps << std::endl;
	}

	// Cleanup
	delete [] lookup;

	return true;
}

template <class T>
bool TomahawkReader::__calculateFST(void){
	this->buffer.resize(this->totempole.getLargestBlockSize() + 1);

	U64 n_read = 0;
		for(U32 i = 0; i < this->totempole.header.blocks; ++i){
			if(!this->nextBlock()){
				std::cerr << "failed to get next block" << std::endl;
				return false;
			}

			// Now have data
			TomahawkIterator<T> controller(this->data.data, this->totempole[i]);
			const Support::TomahawkRunPacked<T>* runs = nullptr;
			const TomahawkEntryMeta<T>* meta = nullptr;

			while(controller.nextVariant(runs, meta)){
				++n_read;
				const U64 hash = XXH64(runs, sizeof(T)*meta->runs, 452930477);
				std::cout << meta->position << '\t' << meta->runs << '\t' << hash << '\n';
			}
		}

		std::cerr << "Total: " << n_read << std::endl;
		return true;
}

template <class T>
bool TomahawkReader::__calculateSFS(void){

}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKREADER_H_ */
