#ifndef TOMAHAWK_TOMAHAWKCALCULATIONS_H_
#define TOMAHAWK_TOMAHAWKCALCULATIONS_H_

#include "TomahawkReader.h"

namespace Tomahawk {

class TomahawkCalculations : public TomahawkReader{
	typedef TomahawkCalculations self_type;

public:
	TomahawkCalculations();
	~TomahawkCalculations();

	bool loadGroups(const std::string& file);
	bool calculateTajimaD(const U32 bin_size);
	bool calculateFST(void);
	bool calculateSFS(void);
	bool calculateIBS(void);
	bool calculateROH(void);
	bool calculateNucleotideDiversity(void);

private:
	template <class T> bool __calculateTajimaD(const U32 bin_size);
	template <class T> bool __calculateFST(void);
	template <class T> bool __calculateSFS(void);
	template <class T> bool __calculateIBS(void);
	template <class T> bool __calculateROH(void);
	template <class T> bool __calculateNucleotideDiversity(void);

private:
	// groups

};

template <class T>
bool TomahawkCalculations::__calculateTajimaD(const U32 bin_size){
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
bool TomahawkCalculations::__calculateFST(void){
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
bool TomahawkCalculations::__calculateSFS(void){
	this->buffer.resize(this->totempole.getLargestBlockSize() + 1);

	// Params
	T* lookup = new T[16];
	std::vector<U64> sfs(2*this->samples, 0);

	for(U32 i = 0; i < this->totempole.header.blocks; ++i){
		if(!this->nextBlock()){
			std::cerr << "failed to get next block" << std::endl;
			return false;
		}

		// Now have data
		TomahawkIterator<T> controller(this->data.data, this->totempole[i]);
		const Support::TomahawkRunPacked<T>* runs = nullptr;
		const TomahawkEntryMeta<T>* meta = nullptr;

		// Cycle over variants in block
		while(controller.nextVariant(runs, meta)){
			// Cycle over runs
			for(U32 i = 0; i < meta->runs; ++i)
				lookup[runs[i].alleles] += runs[i].runs;

			// Frequency for allele 0 and allele 1
			//const U64 f0 = 2*lookup[0] + lookup[1] + lookup[4];
			const U64 f1 = 2*lookup[5] + lookup[1] + lookup[4];

			++sfs[f1];

			// Reset
			memset(lookup, 0, sizeof(T)*16);
		}
	}

	bool folded = true;

	if(!folded){
		std::cout << "AlleleFrequency\tFrequency\n";
		for(U32 i = 0; i < sfs.size(); ++i){
			std::cout << i << '\t' << sfs[i] << '\n';
		}
		std::cout.flush();
	} else {
		std::cout << "AlleleFrequency\tFrequency\n";
		const U32 end = sfs.size() - 1;
		for(U32 i = 0; i < sfs.size()/2; ++i){
			std::cout << i << '\t' << sfs[i] + sfs[end - i] << '\n';
		}
		std::cout.flush();
	}

	delete [] lookup;
	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKCALCULATIONS_H_ */
