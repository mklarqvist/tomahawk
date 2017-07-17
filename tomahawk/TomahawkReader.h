#ifndef TOMAHAWK_TOMAHAWKREADER_H_
#define TOMAHAWK_TOMAHAWKREADER_H_

#include <string>
#include <fstream>
#include <algorithm>
#include <bitset>

#include "MagicConstants.h"
#include "../io/GZController.h"
#include "TomahawkEntryMeta.h"
#include "TotempoleReader.h"
#include "../io/IOConstants.h"
#include "TomahawkCalculateSlave.h"
#include "../interface/Timer.h"
#include "../interface/ProgressBar.h"
#include "TomahawkCalcParameters.h"
#include "../algorithm/Balancer.h"

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
	TomahawkReader(const TotempoleReader& header) :
		samples(0),
		version(0),
		filesize_(0),
		bit_width_(0),
		silent(false),
		threads(std::thread::hardware_concurrency() > 0 ? std::thread::hardware_concurrency() : 1),
		writer(nullptr),
		buffer_(header.getLargestBlockSize()+64),
		data_(this->buffer_),
		outputBuffer_(header.getLargestBlockSize()+64),
		gzip_controller_(),
		totempole_(header),
		balancer()
	{}

	~TomahawkReader(){
		this->buffer_.deleteAll();
		this->data_.deleteAll();
		this->outputBuffer_.deleteAll();
		delete writer;
	}

	bool Open(const std::string input);
	bool ValidateHeader(void);
	bool getBlocks(void);
	bool getBlocks(std::vector<U32>& blocks);
	bool getBlocks(std::vector< std::pair<U32, U32> >& blocks);
	bool getBlock(const U32 blockID);
	bool outputBlocks(std::vector<U32>& blocks);
	bool outputBlocks();

	bool SetR2Threshold(const double min, const double max){
		if(min < 0 || max > 1){
			std::cerr << Helpers::timestamp("ERROR") << "Invalid R-squared values: " << min << '-' << max << std::endl;
			return false;
		}

		if(min > max){
			std::cerr << Helpers::timestamp("ERROR") << "Invalid R-squared values: " << min << '-' << max << std::endl;
			return false;
		}

		this->parameters.R2_min = min - Constants::ALLOWED_ROUNDING_ERROR;
		this->parameters.R2_max = max + Constants::ALLOWED_ROUNDING_ERROR;

		return true;
	}

	bool SetMinimumAlleles(const U64 min){
		this->parameters.minimum_alleles = min;
		return true;
	}

	bool SetThreads(const S32 threads){
		if(threads <= 0){
			std::cerr << Helpers::timestamp("ERROR") << "Invalid number of threads: " << threads << std::endl;
			return false;
		}
		this->threads = threads;
		return true;
	}

	void SetPhased(const bool phased){
		std::cerr << Helpers::timestamp("LOG") << "Forcing phasing of all variants: " << (phased ? "Phased..." : "Unphased...") << std::endl;
		if(phased)
			this->parameters.force = parameter_type::phasedFunction;
		else
			this->parameters.force = parameter_type::unphasedFunction;
	}
	bool SetChunkDesired(const S32 desired){ return(this->balancer.setDesired(desired)); }
	bool SetChunkSelected(const S32 selected){ return(this->balancer.setSelected(selected)); }
	void SetSilent(const bool silent){ this->silent = silent; }

	bool SetPThreshold(const double P){
		if(P > 1 || P < 0){
			std::cerr << Helpers::timestamp("ERROR") << "Invalid P-value threshold: " << P << std::endl;
			return false;
		}
		this->parameters.P_threshold = P;
		return true;
	}

	bool Calculate(std::vector< std::pair<U32,U32> >& blocks);
	bool Calculate(std::vector<U32>& blocks);
	bool Calculate();

	bool SelectWriterOutputType(const IO::TomahawkCalculationWriterInterace::type& writer_type);
	bool OpenWriter(void){
		if(this->writer == nullptr)
			this->SelectWriterOutputType(IO::TomahawkCalculationWriterInterace::type::cout);

		return(this->writer->open());
	}
	void SetOutputType(IO::TomahawkCalculationWriterInterace::compression type){ this->parameters.compression_type = type; }

	bool OpenWriter(const std::string destination){
		if(this->writer == nullptr){
			if(destination == "-"){
				this->SelectWriterOutputType(IO::TomahawkCalculationWriterInterace::type::cout);
				return(this->writer->open());
			}
			this->SelectWriterOutputType(IO::TomahawkCalculationWriterInterace::type::file);
		}

		return(this->writer->open(destination));
	}

private:
	bool __CalculateWrapper();
	template <class T> bool __Calculate();
	void DetermineBitWidth(void);
	U64 GetUncompressedSizes(std::vector<U32>& blocks);
	U64 GetUncompressedSizes(std::vector< std::pair<U32, U32> >& blocks);
	U64 GetUncompressedSizes(void);
	template <class T> bool outputBlock(const U32 blockID);
	template <class T> bool WriteBlock(const char* data, const U32 blockID);
	bool Validate(void);

private:
	U64 samples; // has to match header
	float version; // has to match header
	U64 filesize_;
	BYTE bit_width_;
	std::ifstream stream_; // reader stream
	bool silent;

	//
	U32 threads;
	IO::TomahawkCalculationWriterInterace* writer;
	Interface::ProgressBar progress;

	// Todo: Move out
	IO::BasicBuffer buffer_;
	IO::BasicBuffer data_;
	IO::BasicBuffer outputBuffer_;
	IO::GZController gzip_controller_;
	const TotempoleReader& totempole_;
	Tomahawk::Balancer balancer;

	parameter_type parameters;
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

	U32 readLength = 0;
	if(blockID != this->totempole_.getBlocks() - 1)
		readLength = this->totempole_[blockID + 1].byte_offset - this->totempole_[blockID].byte_offset;
	else
		readLength = this->filesize_ - Constants::eof_length*sizeof(U64) - this->totempole_[this->totempole_.getBlocks()-1].byte_offset;

	if(readLength > this->buffer_.capacity()){
		std::cerr << "impossible: " << readLength << '/' << this->buffer_.capacity() << std::endl;
		exit(1);
	}

	if(!this->stream_.read(&this->buffer_.data[0], readLength)){
		std::cerr << "Failed read: " << this->stream_.good() << '\t' << this->stream_.fail() << '/' << this->stream_.eof() << std::endl;
		std::cerr << this->stream_.gcount() << '/' << readLength << std::endl;
		return false;
	}
	this->buffer_.pointer = readLength;

	char* data_position = &this->data_.data[this->data_.pointer]; // Keep track of position because inflate function moves pointer
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
	TomahawkReaderController<T> tomahawk_controller(data, this->totempole_[blockID]);

	// Foreach variant in block
	for(U32 j = 0; j < tomahawk_controller.support->variants; ++j){
		tomahawk_controller.WriteVariant(this->totempole_, this->outputBuffer_);

		// Next variant
		++tomahawk_controller;

		// Keep flushing regularly
		if(this->outputBuffer_.size() > 65536){
			std::cout.write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
			this->outputBuffer_.reset();
		}
	}

	// Flush last
	std::cout.write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);

	// Reset buffers
	this->outputBuffer_.reset(); // reset
	this->data_.reset(); // reset

	return true;
}

template <class T>
bool TomahawkReader::__Calculate(){
	TomahawkReaderControllerManager<const T> controller(this->totempole_);
	for(U32 i = 0; i < this->blockDataOffsets_.size(); ++i)
		controller.Add(this->blockDataOffsets_[i].data, this->blockDataOffsets_[i].entry);

#if SIMD_AVAILABLE == 1
	std::cerr << Helpers::timestamp("LOG","SIMD") << "Vectorized instructions available: " << SIMD_MAPPING[SIMD_VERSION] << "..." << std::endl;
#else
	std::cerr << Helpers::timestamp("LOG","SIMD") << "No vectorized instructions available..." << std::endl;
#endif
	std::cerr << Helpers::timestamp("LOG","SIMD") << "Building 1-bit representation: ";
	controller.BuildVectorized();
	std::cerr << "Done..." << std::endl;

	const U64 variants = controller.getVariants();
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
	std::cerr << Helpers::timestamp("LOG","CALC") << "Performing " <<  Helpers::ToPrettyString(totalComparisons) << " variant comparisons..."<< std::endl;

	// Setup slaves
	TomahawkCalculateSlave<T>** slaves = new TomahawkCalculateSlave<T>*[this->threads];
	std::vector<std::thread*> thread_pool;

	// Setup workers
	std::cerr << this->parameters << std::endl;
	std::cerr << Helpers::timestamp("LOG","THREAD") << "Spawning " << this->threads << " threads: ";

	for(U32 i = 0; i < this->threads; ++i){
		slaves[i] = new TomahawkCalculateSlave<T>(controller, *this->writer, this->progress, this->parameters, this->balancer.thread_distribution[i]);
		std::cerr << '.';
	}
	std::cerr << std::endl;

	// Start threads
	if(!this->silent)
		this->progress.Start();

	// Todo: fix
	//this->progress.SetDetailed(true);

	// Setup front-end interface
	Interface::Timer timer;
	timer.Start();

	// Begin
	for(U32 i = 0; i < this->threads; ++i)
		thread_pool.push_back(slaves[i]->Start());

	// Join threads
	for(U32 i = 0; i < this->threads; ++i)
		thread_pool[i]->join();

	// Stop progress bar
	this->progress.Stop();

	// Print slave statistics
	std::cerr << Helpers::timestamp("LOG", "THREAD") << "Thread\tPossible\tImpossible\tNoHets\tInsuffucient\tTotal" << std::endl;
	for(U32 i = 0; i < this->threads; ++i)
		std::cerr << Helpers::timestamp("LOG", "THREAD") << i << '\t' << slaves[i]->getPossible() << '\t' << slaves[i]->getImpossible() << '\t' << slaves[i]->getNoHets() << '\t' << slaves[i]->getInsufficientData() << '\t' << slaves[i]->getComparisons() << std::endl;

	// Reduce into first slave
	for(U32 i = 1; i < this->threads; ++i){
		*slaves[0] += *slaves[i];
	}

	std::cerr << Helpers::timestamp("LOG") << "Throughput: " << timer.ElapsedPretty() << " (" << Helpers::ToPrettyString((U64)ceil((double)slaves[0]->getComparisons()/timer.Elapsed().count())) << " pairs of SNP/s, " << Helpers::ToPrettyString((U64)ceil((double)slaves[0]->getComparisons()*this->totempole_.getSamples()/timer.Elapsed().count())) << " genotypes/s)..." << std::endl;
	std::cerr << Helpers::timestamp("LOG") << "Comparisons: " << Helpers::ToPrettyString(slaves[0]->getComparisons()) << " pairwise SNPs and " << Helpers::ToPrettyString(slaves[0]->getComparisons()*this->totempole_.getSamples()) << " pairwise genotypes. Output " << Helpers::ToPrettyString(this->progress.GetOutputCounter()) << "..." << std::endl;

	// Cleanup
	delete [] slaves;

	// Flush writer
	this->writer->flush();

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKREADER_H_ */
