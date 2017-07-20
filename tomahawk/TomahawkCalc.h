#ifndef TOMAHAWK_TOMAHAWKCALC_H_
#define TOMAHAWK_TOMAHAWKCALC_H_

#include "TomahawkReader.h"

namespace Tomahawk {

class TomahawkCalc{
	typedef TomahawkCalc self_type;
	typedef TomahawkCalcParameters parameter_type;
	typedef std::pair<U32,U32> pair_type;
	typedef std::vector<pair_type> pair_vector;
	typedef IO::GenericWriterInterace writer_type;
	typedef Balancer balancer_type;
	typedef TotempoleReader totempole_reader;
	typedef Interface::ProgressBar progress_type;
	typedef TomahawkReader reader_type;

public:
	TomahawkCalc();
	~TomahawkCalc();

	bool Open(const std::string input, const std::string output);

	bool Calculate(pair_vector& blocks);
	bool Calculate(std::vector<U32>& blocks);
	bool Calculate();

	inline parameter_type& getParameters(void){ return(this->parameters); }

private:
	bool OpenWriter(const std::string destination);
	bool SelectWriterOutputType(const writer_type::type writer_type);

	bool CalculateWrapper();
	template <class T> bool Calculate();
	bool WriteTwoHeader(void);
	bool WriteTwoHeaderNatural(void);
	bool WriteTwoHeaderBinary(void);

private:
	bool parameters_validated;
	progress_type progress;
	balancer_type balancer;
	parameter_type parameters;
	reader_type reader;
	writer_type* writer;
};

template <class T>
bool TomahawkCalc::Calculate(){
	const totempole_reader& totempole = this->reader.getTotempole();
	TomahawkBlockManager<const T> controller(totempole);
	for(U32 i = 0; i < totempole.size(); ++i)
		controller.Add(this->reader.getOffsetPair(i).data, this->reader.getOffsetPair(i).entry);

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
	this->progress.SetSamples(totempole.getSamples());

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","CALC") << "Performing " <<  Helpers::ToPrettyString(totalComparisons) << " variant comparisons..."<< std::endl;

	// Setup slaves
	TomahawkCalculateSlave<T>** slaves = new TomahawkCalculateSlave<T>*[this->parameters.n_threads];
	std::vector<std::thread*> thread_pool;

	// Setup workers
	if(!SILENT){
		std::cerr << this->parameters << std::endl;
		std::cerr << Helpers::timestamp("LOG","THREAD") << "Spawning " << this->parameters.n_threads << " threads: ";
	}

	for(U32 i = 0; i < this->parameters.n_threads; ++i){
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
	for(U32 i = 0; i < this->parameters.n_threads; ++i)
		thread_pool.push_back(slaves[i]->Start());

	// Join threads
	for(U32 i = 0; i < this->parameters.n_threads; ++i)
		thread_pool[i]->join();

	// Stop progress bar
	this->progress.Stop();

	// Print slave statistics
	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG", "THREAD") << "Thread\tPossible\tImpossible\tNoHets\tInsuffucient\tTotal" << std::endl;
		for(U32 i = 0; i < this->parameters.n_threads; ++i)
			std::cerr << Helpers::timestamp("LOG", "THREAD") << i << '\t' << slaves[i]->getPossible() << '\t' << slaves[i]->getImpossible() << '\t' << slaves[i]->getNoHets() << '\t' << slaves[i]->getInsufficientData() << '\t' << slaves[i]->getComparisons() << std::endl;
	}

	// Reduce into first slave
	for(U32 i = 1; i < this->parameters.n_threads; ++i)
		*slaves[0] += *slaves[i];

	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG") << "Throughput: " << timer.ElapsedPretty() << " (" << Helpers::ToPrettyString((U64)ceil((double)slaves[0]->getComparisons()/timer.Elapsed().count())) << " pairs of SNP/s, " << Helpers::ToPrettyString((U64)ceil((double)slaves[0]->getComparisons()*totempole.getSamples()/timer.Elapsed().count())) << " genotypes/s)..." << std::endl;
		std::cerr << Helpers::timestamp("LOG") << "Comparisons: " << Helpers::ToPrettyString(slaves[0]->getComparisons()) << " pairwise SNPs and " << Helpers::ToPrettyString(slaves[0]->getComparisons()*totempole.getSamples()) << " pairwise genotypes. Output " << Helpers::ToPrettyString(this->progress.GetOutputCounter()) << "..." << std::endl;
	}

	// Cleanup
	delete [] slaves;

	// Flush writer
	this->writer->flush();

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKCALC_H_ */
