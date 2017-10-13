#ifndef TOMAHAWK_TOMAHAWKCALCULATIONS_H_
#define TOMAHAWK_TOMAHAWKCALCULATIONS_H_

#include "TomahawkReader.h"

namespace Tomahawk {
namespace Stats{
namespace Support{

struct GroupGenotypes{
	GroupGenotypes(void) : count(0){
		memset(&genotypes[0], 0, sizeof(U64)*16);
	}

	~GroupGenotypes(){}

	void reset(void){
		if(count == 0)
			return;

		memset(&genotypes[0], 0, sizeof(U64)*16);
		count = 0;
	}

	void add(const BYTE& genotype, const U64& length){
		this->genotypes[genotype] += length;
		count += length;
	}

	U64 count;
	U64 genotypes[16];
};

template <class T>
struct tajima_helper{
	typedef tajima_helper self_type;

public:
	explicit tajima_helper(void) :
		n(0),
		a1(0), a2(0), b1(0), b2(0), c1(0), c2(0), e1(0), e2(0),
		s_minor(0),
		s_major(0),
		k_hat(0),
		pi_cum(0),
		f0(0),
		f1(0),
		n_snps(0),
		D(0),
		meanMAF(0),
		meanPI(0)
	{
	}

	void operator()(const double n){
		this->n = 2*n;

		for(U32 i = 1; i < this->n - 1; ++i){
			this->a1 += (double)1/i;
			this->a2 += (double)1/(i*i);
		}

		this->b1 = double(this->n+1) / 3.0 / double(this->n-1);
		this->b2 = 2.0 * double(this->n*this->n + this->n + 3) / 9.0 / double(this->n) / double(this->n-1);
		this->c1 = this->b1 - (1.0 / this->a1);
		this->c2 = this->b2 - (double(this->n+2)/double(this->a1*this->n)) + (this->a2/this->a1/this->a1);
		this->e1 = this->c1 / this->a1;
		this->e2 = this->c2 / ((this->a1*this->a1) + this->a2);
	}

	~tajima_helper(){}

	void reset(void){
		this->n_snps = 0;
		this->k_hat = 0; this->pi_cum = 0;
		this->s_minor = 0; this->s_major = 0;
	}

	void operator++(void){ ++this->n_snps; }

	void update(const GroupGenotypes& g){
		// Frequency for allele 0 and allele 1
		this->f0 = 2*g.genotypes[0] + g.genotypes[1] + g.genotypes[4];
		this->f1 = 2*g.genotypes[5] + g.genotypes[1] + g.genotypes[4];
		const U32 total = this->f0 + this->f1;

		if(total == 0)
			return;

		// Update minor allele frequency
		if(this->f0 < this->f1){
			this->s_minor += (double)this->f0/total;
			this->s_major += (double)this->f1/total;
		} else {
			this->s_minor += (double)this->f1/total;
			this->s_major += (double)this->f0/total;
		}

		// Update k_hat
		this->k_hat += (double)this->f0/total*(double)this->f1/total;
		++this->n_snps;

		// Update nucleotide diversity
		this->pi_cum += 2*this->f0*this->f1 / ((double)total * (total - 1));
	}

	void calc(void){
		if(this->n_snps == 0){
			this->D = 0;
			this->meanMAF = 0;
			this->meanPI = 0;
			return;
		}

		const double S = this->n_snps;
		const double pi = 2.0*this->k_hat*this->n/double(this->n-1);
		const double tw = double(S) / this->a1;
		const double var = (this->e1*S) + this->e2*S*(S-1);
		this->D = (pi - tw) / sqrt(var);
		this->meanMAF = ((double)this->s_minor / (this->s_minor + this->s_major)) / S;
		this->meanPI = pi_cum/S;
	}

public:
	double n;
	double a1, a2, b1, b2, c1, c2, e1, e2;
	double s_minor, s_major;
	double k_hat;
	double pi_cum;
	U64 f0, f1;
	U64 n_snps;

	// Output
	double D, meanMAF, meanPI;
};

}
}

class TomahawkCalculations : public TomahawkReader{
	typedef TomahawkCalculations self_type;
	typedef Tomahawk::Hash::HashTable<std::string, S32> hash_table;
	typedef std::vector<U64> occ_vector;
	typedef std::vector<occ_vector> occ_matrix;

public:
	struct GroupPair{
		GroupPair(const std::string& name) : n_entries(1), name(name){}
		~GroupPair(){}

		void operator++(void){ ++this->n_entries; }
		void operator--(void){ --this->n_entries; }
		void operator+=(const U32 p){ this->n_entries += p; }
		void operator-=(const U32 p){ this->n_entries -= p; }

		U32 n_entries;
		std::string name;
	};

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
	template <class T> bool __calculateTajimaDGrouped(const U32 bin_size);
	template <class T> bool __calculateFST(void);
	template <class T> bool __calculateSFS(void);
	template <class T> bool __calculateSFSGrouped(void);
	template <class T> bool __calculateIBS(void);
	template <class T> bool __calculateROH(void);
	template <class T> bool __calculateNucleotideDiversity(void);

private:
	// groups
	occ_matrix Occ;
	std::vector<GroupPair> groups;
	hash_table* group_htable;
};

template <class T>
bool TomahawkCalculations::__calculateTajimaDGrouped(const U32 bin_size){
	if(this->Occ.size() == 0){
		std::cerr << "no occ matrix loaded" << std::endl;
		return false;
	}

	this->buffer.resize(this->totempole.getLargestBlockSize() + 1);

	// Params
	std::vector<Stats::Support::GroupGenotypes> genotypes(this->Occ[0].size());
	std::vector< Stats::Support::tajima_helper<T> > helpers(genotypes.size());
	for(U32 i = 0; i < genotypes.size(); ++i)
		helpers[i](this->groups[i].n_entries);

	U64 current_bin = 0;
	U32 prevPos = 0;
	S32 previous_contigID = 0;
	U64 cumSum = 0;

	std::cout << "n_snps\tcontigID\tbinFrom\tbinTo\tcumBinFrom\tcumBinTo\tminD\tmeanD\tmaxD\t";
	for(U32 i = 0; i < genotypes.size() - 1; ++i)
		std::cout << this->groups[i].name << '\t';
	std::cout << this->groups[genotypes.size() - 1].name << std::endl;

	for(U32 i = 0; i < this->totempole.header.blocks; ++i){
		if(!this->nextBlock()){
			std::cerr << "failed to get next block" << std::endl;
			return false;
		}

		// Reset arrays
		for(U32 i = 0; i < genotypes.size(); ++i)
			genotypes[i].reset();

		// Now have data
		TomahawkIterator<T> controller(this->data.data, this->totempole[i]);
		const Support::TomahawkRunPacked<T>* runs = nullptr;
		const TomahawkEntryMeta<T>* meta = nullptr;

		// Cycle over variants
		while(controller.nextVariant(runs, meta)){
			// Assert file is sorted
			// This should always be true
			if((meta->position < prevPos && previous_contigID == this->totempole[i].contigID) || previous_contigID > this->totempole[i].contigID){
				std::cerr << "unsorted" << std::endl;
				exit(1);
			}

			// Count number of genotypes / group
			U64 cumPos = 0;
			for(U32 i = 0; i < meta->runs; ++i){
				for(U32 k = 0; k < this->Occ[0].size(); ++k)
					genotypes[k].add(runs[i].alleles, this->Occ[cumPos + runs[i].runs][k] - this->Occ[cumPos][k]);

				cumPos += runs[i].runs;
			}

			// If we reach the end of a bin or switch chromosome
			// then output calculations
			const U32 bin = (U32)((double)meta->position/bin_size)*bin_size;
			if(bin != current_bin || this->totempole[i].contigID != previous_contigID){
				// calc
				// do something with the output
				double minD = 1000000;
				double maxD = -1000000;
				double sumD = 0;
				U64 n_D = 0;
				U32 toPos = current_bin + bin_size;
				if(toPos > this->totempole.contigs[previous_contigID].maxPosition)
					toPos = this->totempole.contigs[previous_contigID].maxPosition;

				for(U32 i = 0; i < genotypes.size(); ++i){
					helpers[i].calc();
					if(helpers[i].n_snps == 0)
						continue;

					sumD += helpers[i].D;
					if(helpers[i].D < minD) minD = helpers[i].D;
					if(helpers[i].D > maxD) maxD = helpers[i].D;
					if(helpers[i].n_snps > 0) ++n_D;
				}

				if(n_D > 0){
				std::cout << n_D << '\t' << previous_contigID << '\t' << current_bin << '\t' << toPos << '\t'
						<< cumSum << '\t' << cumSum + (toPos - current_bin) << '\t' << minD << '\t'
						<< sumD/n_D << '\t' << maxD << '\t';
				}

				if(n_D > 0){
					for(U32 i = 0; i < genotypes.size() - 1; ++i){
						std::cout << helpers[i].D << '\t';
						helpers[i].reset();
					}
					std::cout << helpers[genotypes.size() - 1].D << '\n';

				} else {
					// Cleanup
					for(U32 i = 0; i < genotypes.size(); ++i)
						helpers[i].reset();
				}

				// Update cumsum (used for plotting)
				cumSum += toPos - current_bin;

				// Reset
				previous_contigID = this->totempole[i].contigID;
			}

			// Update data
			for(U32 i = 0; i < genotypes.size(); ++i)
				helpers[i].update(genotypes[i]);

			// Updates
			prevPos = meta->position;
			current_bin = (U32)((double)meta->position/bin_size)*bin_size;

			// Reset arrays
			for(U32 i = 0; i < genotypes.size(); ++i)
				genotypes[i].reset();
		}
	}

	return true;
}

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
	T lookup[16];
	U32 counts = 0;

	U64 f0 = 0;
	U64 f1 = 0;
	double k_hat = 0;
	double pi_cum = 0;

	U32 current_bin = 0;
	U32 prevPos = 0;
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
			const U32 bin = (U32)((double)meta->position/bin_size)*bin_size;
			if(bin != current_bin || this->totempole[i].contigID != previous_contigID){
				if(counts > 0){
					U32 toPos = current_bin + bin_size;
					if(toPos > this->totempole.contigs[previous_contigID].maxPosition)
						toPos = this->totempole.contigs[previous_contigID].maxPosition;

					const double S = counts;
					const double pi = 2.0*k_hat*n/double(n-1);
					const double tw = double(S) / a1;
					const double var = (e1*S) + e2*S*(S-1);
					const double D = (pi - tw) / sqrt(var);
					const double meanMAF = ((double)s_minor / (s_minor + s_major)) / S;

					std::cout << S << '\t' << previous_contigID << '\t' << current_bin << '\t' << toPos << '\t'
							  << cumSum << '\t' << cumSum + (toPos - current_bin) << '\t' << D << '\t' << pi_cum/S << '\t' << meanMAF << '\t' << k_hat << std::endl;
					cumSum += bin_size;
				}

				// Reset
				counts = 0;
				f0 = 0; f1 = 0;
				k_hat = 0; pi_cum = 0;
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
			current_bin = (U32)((double)meta->position/bin_size)*bin_size;

			// Reset
			memset(lookup, 0, sizeof(T)*16);
		}
	}

	// Output last partition
	if(counts > 0){
		U32 toPos = current_bin + bin_size;
		if(toPos > this->totempole.contigs[previous_contigID].maxPosition)
			toPos = this->totempole.contigs[previous_contigID].maxPosition;

		const double S = counts;
		const double pi = 2.0*k_hat*n/double(n-1);
		const double tw = double(S) / a1;
		const double var = (e1*S) + e2*S*(S-1);
		const double D = (pi - tw) / sqrt(var);
		const double meanMAF = ((double)s_minor / (s_minor + s_major)) / S;

		std::cout << S << '\t' << previous_contigID << '\t' << current_bin << '\t' << toPos << '\t'
				  << cumSum << '\t' << cumSum + (toPos - current_bin) << '\t' << D << '\t' << pi_cum/S << '\t' << meanMAF << '\t' << k_hat << std::endl;
		cumSum += bin_size;
	}
	return true;
}

template <class T>
bool TomahawkCalculations::__calculateFST(void){
	this->buffer.resize(this->totempole.getLargestBlockSize() + 1);

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
			//
		}
	}

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

template <class T>
bool TomahawkCalculations::__calculateSFSGrouped(void){
	if(this->Occ.size() == 0){
		std::cerr << "no occ matrix loaded" << std::endl;
		return false;
	}

	this->buffer.resize(this->totempole.getLargestBlockSize() + 1);

	// Params
	std::vector<Stats::Support::GroupGenotypes> genotypes(this->Occ[0].size());
	std::vector< std::vector<U64> > sfs(this->Occ[0].size(), std::vector<U64>(2*this->samples,0));
	U64 cumPos = 0;

	// Get reference group name
	const std::string ref = "EUR"; // Todo: make external
	S32* ref_groupID = nullptr;
	if(!this->group_htable->GetItem(&ref[0], &ref, ref_groupID, ref.length())){
		std::cerr << "ref: " << ref << " not defined" << std::endl;
		return false;
	}
	std::cerr << "ref is: " << *ref_groupID << std::endl;
	const S32 refID = *ref_groupID;
	const U32 expected_n_ref = this->groups[refID].n_entries;

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
			for(U32 i = 0; i < meta->runs; ++i){
				for(U32 k = 0; k < this->Occ[0].size(); ++k)
					genotypes[k].add(runs[i].alleles, this->Occ[cumPos + runs[i].runs][k] - this->Occ[cumPos][k]);

				cumPos += runs[i].runs;
			}

			const U32 ref_counts = genotypes[refID].genotypes[0] + genotypes[refID].genotypes[1] + genotypes[refID].genotypes[4] + genotypes[refID].genotypes[5];
			if(ref_counts == 0){
				cumPos = 0;

				// Reset
				for(U32 k = 0; k < this->Occ[0].size(); ++k)
					genotypes[k].reset();

				continue;
			}

			if(genotypes[refID].genotypes[1] > 0 || genotypes[refID].genotypes[4] > 0)
				std::cerr << genotypes[refID].genotypes[0] << '\t' << genotypes[refID].genotypes[1] << '\t' << genotypes[refID].genotypes[4] << '\t' << genotypes[refID].genotypes[5] << std::endl;

			// Reset and update
			for(U32 k = 0; k < this->Occ[0].size(); ++k){
				const U64 f0 = 2*genotypes[k].genotypes[0] + genotypes[k].genotypes[1] + genotypes[k].genotypes[4];
				const U64 f1 = 2*genotypes[k].genotypes[5] + genotypes[k].genotypes[1] + genotypes[k].genotypes[4];

				if(genotypes[refID].genotypes[0] == expected_n_ref)
					++sfs[k][f1];
				else if(genotypes[refID].genotypes[5] == expected_n_ref)
					++sfs[k][f0];
				else
					continue;
			}

			// Reset
			for(U32 k = 0; k < this->Occ[0].size(); ++k)
				genotypes[k].reset();

			cumPos = 0;
		}
	}

	// Output
	// Header
	for(U32 i = 0; i < this->groups.size(); ++i){
		std::cout << this->groups[i].name << '\t';
	}
	std::cout << std::endl;

	for(U32 i = 0; i < sfs[0].size(); ++i){
		for(U32 k = 0; k < this->Occ[0].size(); ++k){
			std::cout << sfs[k][i] << '\t';
		}
		std::cout << std::endl;
	}

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKCALCULATIONS_H_ */
