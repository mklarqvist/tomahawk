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

	// Used to keep track of char pointer offsets in buffer
	// and what totempole entry is associated with that position
	struct DataOffsetPair{
		DataOffsetPair(const char* data, const TotempoleEntry& entry) : entry(entry), data(data){}
		~DataOffsetPair(){}

		const TotempoleEntry& entry;
		const char* data;
	};

public:
	TomahawkReader(const TotempoleReader& header);
	~TomahawkReader();

	// Reader functions
	bool Open(const std::string input);
	bool ValidateHeader(void);
	bool getBlocks(void);
	bool getBlocks(std::vector<U32>& blocks);
	bool getBlocks(std::vector< std::pair<U32, U32> >& blocks);
	bool getBlock(const U32 blockID);
	bool outputBlocks(std::vector<U32>& blocks);
	bool outputBlocks();

	// Calc functions
	bool SetR2Threshold(const double min, const double max);
	bool SetMinimumAlleles(const U64 min);
	bool SetThreads(const S32 threads);
	void SetPhased(const bool phased);
	bool SetChunkDesired(const S32 desired){ return(this->balancer.setDesired(desired)); }
	bool SetChunkSelected(const S32 selected){ return(this->balancer.setSelected(selected)); }
	bool SetPThreshold(const double P);
	void setDetailedProgress(const bool yes);

	bool Calculate(std::vector< std::pair<U32,U32> >& blocks);
	bool Calculate(std::vector<U32>& blocks);
	bool Calculate();
	bool SelectWriterOutputType(const IO::GenericWriterInterace::type writer_type);
	void SetOutputType(IO::GenericWriterInterace::compression type){ this->parameters.compression_type = type; }
	bool OpenWriter(void);
	bool OpenWriter(const std::string destination);
	inline const BYTE& getBitWidth(void) const{ return(this->bit_width_); }
	inline const TotempoleReader& getTotempole(void) const{ return(this->totempole_); }

private:
	bool WriteTwoHeader(void);
	bool WriteTwoHeaderNatural(void);
	bool WriteTwoHeaderBinary(void);

private:
	bool __CalculateWrapper();
	void DetermineBitWidth(void);
	U64 GetUncompressedSizes(std::vector<U32>& blocks);
	U64 GetUncompressedSizes(std::vector< std::pair<U32, U32> >& blocks);
	U64 GetUncompressedSizes(void);

	// Calc functions
	template <class T> bool __Calculate();
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
	//bool silent;

	//
	IO::GenericWriterInterace* writer;

	U32 threads;
	Interface::ProgressBar progress;
	const TotempoleReader& totempole_;
	Tomahawk::Balancer balancer;
	parameter_type parameters;

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

template <class T>
bool TomahawkReader::__Calculate(){
	TomahawkBlockManager<const T> controller(this->totempole_);
	for(U32 i = 0; i < this->blockDataOffsets_.size(); ++i)
		controller.Add(this->blockDataOffsets_[i].data, this->blockDataOffsets_[i].entry);

	if(!SILENT){
#if SIMD_AVAILABLE == 1
	std::cerr << Helpers::timestamp("LOG","SIMD") << "Vectorized instructions available: " << SIMD_MAPPING[SIMD_VERSION] << "..." << std::endl;
#else
	std::cerr << Helpers::timestamp("LOG","SIMD") << "No vectorized instructions available..." << std::endl;
#endif
	std::cerr << Helpers::timestamp("LOG","SIMD") << "Building 1-bit representation: ";
	}

	// Build 1-bit representation
	controller.BuildVectorized();

	if(!SILENT)
		std::cerr << "Done..." << std::endl;

	const U64 variants = controller.getVariants();

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","CALC") << "Total " << Helpers::ToPrettyString(variants) << " variants..." << std::endl;

	// Todo: validate
	U64 totalComparisons = 0;
	for(U32 i = 0; i < this->balancer.thread_distribution.size(); ++i){
		for(U32 j = 0; j < this->balancer.thread_distribution[i].size(); ++j){
			//std::cerr << this->balancer.thread_distribution[i][j] << ':' << std::endl;
			if(this->balancer.thread_distribution[i][j].staggered){
				for(U32 from = this->balancer.thread_distribution[i][j].fromRow; from < this->balancer.thread_distribution[i][j].toRow; ++from){
					for(U32 col = from; col < this->balancer.thread_distribution[i][j].toColumn; ++col){
						//std::cerr << '\t' << from << ":" << col << '\t';
						if(from == col){
							const U32 size = controller[from].size();
							totalComparisons += (size*size - size)/2;
							//std::cerr << (size*size - size)/2 << std::endl;
						} else {
							totalComparisons += controller[from].size() * controller[col].size();
							//std::cerr << controller[from].size() * controller[col].size() << std::endl;
						}
					}
				}
			} else {
				for(U32 from = this->balancer.thread_distribution[i][j].fromRow; from < this->balancer.thread_distribution[i][j].toRow; ++from){
					for(U32 col = this->balancer.thread_distribution[i][j].fromColumn; col < this->balancer.thread_distribution[i][j].toColumn; ++col){
						//std::cerr << '\t' << from << ":" << col << '\t';
						if(from == col){
							const U32 size = controller[from].size();
							totalComparisons += (size*size - size)/2;
							//std::cerr << (size*size - size)/2 << std::endl;
						} else {
							totalComparisons += controller[from].size() * controller[col].size();
							//std::cerr << controller[from].size() * controller[col].size() << std::endl;
						}
					}
				}
			}
		}
	}
	this->progress.SetComparisons(totalComparisons);

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","CALC") << "Performing " <<  Helpers::ToPrettyString(totalComparisons) << " variant comparisons..."<< std::endl;

	// Setup slaves
	TomahawkCalculateSlave<T>** slaves = new TomahawkCalculateSlave<T>*[this->threads];
	std::vector<std::thread*> thread_pool;

	// Setup workers
	if(!SILENT){
		std::cerr << this->parameters << std::endl;
		std::cerr << Helpers::timestamp("LOG","THREAD") << "Spawning " << this->threads << " threads: ";
	}

	for(U32 i = 0; i < this->threads; ++i){
		slaves[i] = new TomahawkCalculateSlave<T>(controller, *this->writer, this->progress, this->parameters, this->balancer.thread_distribution[i]);
		if(!SILENT)
			std::cerr << '.';
	}
	if(!SILENT)
		std::cerr << std::endl;

	// Start threads
	if(!SILENT)
		this->progress.Start();

	// Setup front-end interface
	Interface::Timer timer;
	timer.Start();

	// Write TWO output header
	this->WriteTwoHeader();

	// Begin
	for(U32 i = 0; i < this->threads; ++i)
		thread_pool.push_back(slaves[i]->Start());

	// Join threads
	for(U32 i = 0; i < this->threads; ++i)
		thread_pool[i]->join();

	// Stop progress bar
	this->progress.Stop();

	// Print slave statistics
	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG", "THREAD") << "Thread\tPossible\tImpossible\tNoHets\tInsuffucient\tTotal" << std::endl;
		for(U32 i = 0; i < this->threads; ++i)
			std::cerr << Helpers::timestamp("LOG", "THREAD") << i << '\t' << slaves[i]->getPossible() << '\t' << slaves[i]->getImpossible() << '\t' << slaves[i]->getNoHets() << '\t' << slaves[i]->getInsufficientData() << '\t' << slaves[i]->getComparisons() << std::endl;
	}

	// Reduce into first slave
	for(U32 i = 1; i < this->threads; ++i)
		*slaves[0] += *slaves[i];

	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG") << "Throughput: " << timer.ElapsedPretty() << " (" << Helpers::ToPrettyString((U64)ceil((double)slaves[0]->getComparisons()/timer.Elapsed().count())) << " pairs of SNP/s, " << Helpers::ToPrettyString((U64)ceil((double)slaves[0]->getComparisons()*this->totempole_.getSamples()/timer.Elapsed().count())) << " genotypes/s)..." << std::endl;
		std::cerr << Helpers::timestamp("LOG") << "Comparisons: " << Helpers::ToPrettyString(slaves[0]->getComparisons()) << " pairwise SNPs and " << Helpers::ToPrettyString(slaves[0]->getComparisons()*this->totempole_.getSamples()) << " pairwise genotypes. Output " << Helpers::ToPrettyString(this->progress.GetOutputCounter()) << "..." << std::endl;
	}

	// Cleanup
	delete [] slaves;

	// Flush writer
	this->writer->flush();

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKREADER_H_ */
