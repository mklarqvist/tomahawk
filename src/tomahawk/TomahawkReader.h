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
	bool loadGroups(const std::string& file);
	bool calculateTajimaD(const U32 bin_size);
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
	template <class T> bool __calculateTajimaD(const U32 bin_size);
	template <class T> bool __calculateFST(void);
	template <class T> bool __calculateSFS(void);

private:
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

template <class T>
bool TomahawkReader::__calculateTajimaD(const U32 bin_size){
	this->buffer.resize(this->totempole.getLargestBlockSize() + 1);

	// Constants
	const double n = 2*this->totempole.header.samples;
	double a1 = 0, a2 = 0;
	for(U32 i = 1; i < n - 1; ++i){
		a1 += (double)1/i;
		a2 += (double)1/(i*i);
	}

	double b1 = double(n+1) / 3.0 / double(n-1);
	double b2 = 2.0 * double(n*n + n + 3) / 9.0 / double(n) / double(n-1);
	double c1 = b1 - (1.0 / a1);
	double c2 = b2 - (double(n+2)/double(a1*n)) + (a2/a1/a1);
	double e1 = c1 / a1;
	double e2 = c2 / ((a1*a1) + a2);

	// Variables
	T* lookup = new T[16];
	U32 counts = 0;
	U32 prevPos = 0;
	U64 f0 = 0;
	U64 f1 = 0;
	double k_hat = 0;
	double pi_cum = 0;

	U32 current_bin = 0;
	S32 previous_contigID = 0;
	U64 cumSum = 0;
	double s_minor = 0, s_major = 0;

	std::cout << "n_snps\tcontigID\tbinFrom\tbinTo\tcumBinFrom\tcumBinTo\tTajimaD\tmeanPI\tmeanMAF\tk_hat" << std::endl;
	for(U32 i = 0; i < this->totempole.header.blocks; ++i){
		if(!this->nextBlock()){
			std::cerr << "failed to get next block" << std::endl;
			return false;
		}

		// Reset array
		memset(lookup, 0, sizeof(T)*16);

		// Now have data
		TomahawkIterator<T> controller(this->data.data, this->totempole[i]);
		const Support::TomahawkRunPacked<T>* runs = nullptr;
		const TomahawkEntryMeta<T>* meta = nullptr;

		while(controller.nextVariant(runs, meta)){
			// Assert file is sorted
			// This should always be true
			if((meta->position < prevPos && previous_contigID == this->totempole[i].contigID) || previous_contigID > this->totempole[i].contigID){
				std::cerr << "unsorted" << std::endl;
				exit(1);
			}

			// Count number of genotypes
			for(U32 i = 0; i < meta->runs; ++i)
				lookup[runs[i].alleles] += runs[i].runs;

			// If we reach the end of a bin or switch chromosome
			// then output calculations
			if((meta->position/bin_size)*bin_size != current_bin || this->totempole[i].contigID != previous_contigID){
				if(counts > 0){
					const double S = counts;
					const double pi = 2.0*k_hat*n/double(n-1);
					const double tw = double(S) / a1;
					const double var = (e1*S) + e2*S*(S-1);
					const double D = (pi - tw) / sqrt(var);
					const double meanMAF = ((double)s_minor / (s_minor + s_major)) / S;

					std::cout << counts << '\t' << previous_contigID << '\t' << current_bin << '\t' << current_bin + bin_size << '\t' << cumSum << '\t' << cumSum+bin_size << '\t' << D << '\t' << pi_cum/S << '\t' << meanMAF << '\t' << k_hat << std::endl;
					cumSum += bin_size;
				}

				// Reset
				counts = 0;
				f0 = 0; f1 = 0;
				k_hat = 0; pi_cum = 0;
				current_bin = (meta->position/bin_size)*bin_size;
				previous_contigID = this->totempole[i].contigID;
				s_minor = 0; s_major = 0;
			}

			// Frequency for allele 0 and allele 1
			f0 = 2*lookup[0] + lookup[1] + lookup[4];
			f1 = 2*lookup[5] + lookup[1] + lookup[4];
			const U32 total = f0 + f1;

			// Update minor allele frequency
			if(f0 < f1){ s_minor += (double)f0/total; s_major += (double)f1/total; }
			else       { s_minor += (double)f1/total; s_major += (double)f0/total; }

			// Update k_hat
			k_hat += (double)f0/total*(double)f1/total;
			++counts;

			// Update nucleotide diversity
			pi_cum += 2*f0*f1 / ((double)total * (total - 1));

			// Updates
			prevPos = meta->position;

			// Reset
			memset(lookup, 0, sizeof(T)*16);
		}
	}

	// Output last partition
	if(counts > 0){
		const double S = counts;
		const double pi = 2.0*k_hat*n/double(n-1);
		const double tw = double(S) / a1;
		const double var = (e1*S) + e2*S*(S-1);
		const double D = (pi - tw) / sqrt(var);
		const double meanMAF = ((double)s_minor / (s_minor + s_major)) / S;

		std::cout << counts << '\t' << previous_contigID << '\t' << current_bin << '\t' << current_bin + bin_size << '\t' << cumSum << '\t' << cumSum+bin_size << '\t' << D << '\t' << pi_cum/S << '\t' << meanMAF << '\t' << k_hat << std::endl;
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
