#ifndef TOMAHAWK_TOMAHAWKCALC_H_
#define TOMAHAWK_TOMAHAWKCALC_H_

#include "TomahawkReader.h"
#include "twk_reader_implementation.h"
#include "genotype_meta_container_reference.h"
#include "../io/output_writer.h"
#include "../index/index.h"
#include "../algorithm/load_balancer_ld.h"

namespace Tomahawk {

class TomahawkCalc{
	typedef TomahawkCalc               self_type;
	typedef TomahawkCalcParameters     parameter_type;
	typedef std::pair<U32,U32>         pair_type;
	typedef std::vector<pair_type>     pair_vector;
	typedef LoadBalancerLD             balancer_type;
	typedef TomahawkHeader             header_type;
	typedef Interface::ProgressBar     progress_type;
	typedef TomahawkReader             reader_type;

public:
	TomahawkCalc();
	~TomahawkCalc();
	bool Open(const std::string input, const std::string output);
	bool Calculate(pair_vector& blocks);
	bool Calculate(std::vector<U32>& blocks);
	bool Calculate();
	inline parameter_type& getParameters(void){ return(this->parameters); }

private:
	template <class T> bool Calculate();

private:
	std::string    input_file;
	std::string    output_file;
	bool           parameters_validated;
	progress_type  progress;
	balancer_type  balancer;
	parameter_type parameters;
	reader_type    reader;
};

template <class T>
bool TomahawkCalc::Calculate(){
	// Retrieve reference to Totempole reader
	header_type& header = this->reader.getHeader();
	header.addLiteral("\n##tomahawk_calcCommand=" + Helpers::program_string());
	header.addLiteral("\n##tomahawk_calcInterpretedCommand=" + this->parameters.getInterpretedString());


	IO::OutputWriter writer;
	if(!writer.open(this->output_file)){
		std::cerr << Helpers::timestamp("ERROR", "TWI") << "Failed to open..." << std::endl;
		return false;
	}

	writer.writeHeaders(this->reader.getHeader());


	if(!SILENT){
	#if SIMD_AVAILABLE == 1
		std::cerr << Helpers::timestamp("LOG","SIMD") << "Vectorized instructions available: " << SIMD_MAPPING[SIMD_VERSION] << "..." << std::endl;
	#else
		std::cerr << Helpers::timestamp("LOG","SIMD") << "No vectorized instructions available..." << std::endl;
	#endif
		std::cerr << Helpers::timestamp("LOG","SIMD") << "Building 1-bit representation: ";
	}

	// Construct Tomahawk manager
	//TomahawkReaderImpl<T> impl(totempole.header.samples, this->reader.DataOffsetSize()+1);

	//std::cerr << "not implemented" << std::endl;
	//exit(1);

	GenotypeMetaContainerReference<T> references(header.magic_.getNumberSamples(), this->reader.DataOffsetSize()+1);

	U64 n_variants = 0;
	for(U32 i = 0; i < this->reader.DataOffsetSize(); ++i){
		// Hard copy of data into STL-like containers
		/*
		impl.addDataBlock(this->reader.getOffsetPair(i).data,
                          this->reader.getOffsetPair(i).l_buffer,
                          this->reader.getOffsetPair(i).entry);
		*/

		// Reference interpretation of char buffer in psuedo-iterator containers
		// directly from unaligned memory
		references.addDataBlock(this->reader.getOffsetPair(i).data,
                                this->reader.getOffsetPair(i).l_buffer,
                                this->reader.getOffsetPair(i).entry);

		n_variants += references[i].getTotempole().size();
	}

	if(!SILENT)
		std::cerr << "Done..." << std::endl;

	// Number of variants in memory
	const U64 variants = references.countVariants();

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","CALC") << "Total " << Helpers::ToPrettyString(variants) << " variants..." << std::endl;

	// Update progress bar with data
	this->progress.SetComparisons(this->balancer.n_comparisons_chunk);
	this->progress.SetSamples(header.magic_.getNumberSamples());
	this->progress.SetDetailed(this->parameters.detailed_progress);

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","CALC") << "Performing " <<  Helpers::ToPrettyString(this->balancer.n_comparisons_chunk) << " variant comparisons..."<< std::endl;

	// Setup slaves
	LDSlave<T>** slaves = new LDSlave<T>*[this->parameters.n_threads];
	std::vector<std::thread*> thread_pool;

	if(!SILENT){
		if(this->parameters.fast_mode){
			std::cerr << Helpers::timestamp("LOG") << "Running in fast mode! No matrices will be built..." << std::endl;
		} else {
			std::cerr << Helpers::timestamp("LOG") << "Running in complete mode! 2x2/3x3/4x4 matrices will be built..." << std::endl;
		}

	}

	// Setup workers
	if(!SILENT){
		std::cerr << this->parameters << std::endl;
		std::cerr << Helpers::timestamp("LOG","THREAD") << "Spawning " << this->parameters.n_threads << " threads: ";
	}

	for(U32 i = 0; i < this->parameters.n_threads; ++i){
		slaves[i] = new LDSlave<T>(references, writer, this->progress, this->parameters, this->balancer.thread_distribution[i]);
		if(!SILENT)
			std::cerr << '.';
	}

	if(!SILENT)
		std::cerr << std::endl;

	// Start progress tracker
	if(!SILENT)
		this->progress.Start();

	// Setup front-end interface
	Interface::Timer timer;
	timer.Start();

	// Start workers
	for(U32 i = 0; i < this->parameters.n_threads; ++i)
		thread_pool.push_back(slaves[i]->Start());

	// Join workers
	for(U32 i = 0; i < this->parameters.n_threads; ++i)
		thread_pool[i]->join();

	// Stop progress tracker
	this->progress.Stop();

	// Print slave statistics
	/*
	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG", "THREAD") << "Thread\tPossible\tImpossible\tNoHets\tInsuffucient\tTotal" << std::endl;
		for(U32 i = 0; i < this->parameters.n_threads; ++i)
			std::cerr << Helpers::timestamp("LOG", "THREAD") << i << '\t' << slaves[i]->getPossible() << '\t' << slaves[i]->getImpossible() << '\t' << slaves[i]->getNoHets() << '\t' << slaves[i]->getInsufficientData() << '\t' << slaves[i]->getComparisons() << std::endl;
	}
	*/

	// Reduce into first slave
	for(U32 i = 1; i < this->parameters.n_threads; ++i)
		*slaves[0] += *slaves[i];

	writer = slaves[0]->getWriter();

	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG") << "Throughput: " << timer.ElapsedString() << " (" << Helpers::ToPrettyString((U64)ceil((double)this->balancer.n_comparisons_chunk/timer.Elapsed().count())) << " pairs of SNP/s, " << Helpers::ToPrettyString((U64)ceil((double)this->balancer.n_comparisons_chunk*header.magic_.getNumberSamples()/timer.Elapsed().count())) << " genotypes/s)..." << std::endl;
		std::cerr << Helpers::timestamp("LOG") << "Comparisons: " << Helpers::ToPrettyString(this->balancer.n_comparisons_chunk) << " pairwise SNPs and " << Helpers::ToPrettyString(this->balancer.n_comparisons_chunk*header.magic_.getNumberSamples()) << " pairwise genotypes..." << std::endl;
		std::cerr << Helpers::timestamp("LOG") << "Output: " << Helpers::ToPrettyString(writer.sizeEntries()) << " entries into " << Helpers::ToPrettyString(writer.sizeBlocks()) << " blocks..." << std::endl;
	}

	// Cleanup
	delete [] slaves;

	// Flush writer
	writer.flush();
	writer.writeFinal();

	/*
	if(!writer.finalise()){
		std::cerr << Helpers::timestamp("ERROR", "INDEX") << "Failed to finalize..." << std::endl;
		return false;
	}

	writer.close();
	*/

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKCALC_H_ */
