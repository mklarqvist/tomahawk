#ifndef TOMAHAWKIMPORTRLE_H_
#define TOMAHAWKIMPORTRLE_H_

#include <algorithm>
#include <bitset>

#include "io/bcf/bcf_reader.h"
#include "io/vcf/vcf_lines.h"
#include "math/fisher_math.h"
#include "tomahawk/meta_entry.h"
#include "tomahawk/import_filters.h"

namespace tomahawk{
namespace algorithm{

// Todo: These should all be moved to VCF or BCF entry class
struct TomahawkImportRLEHelper{
	typedef TomahawkImportRLEHelper self_type;

	TomahawkImportRLEHelper(const U64 expectedSamples) :
		MAF(0),
		MGF(0),
		HWE_P(0),
		missingValues(0),
		phased(false),
		expectedSamples(expectedSamples)
	{
		memset(&this->countsGenotypes[0], 0, sizeof(U64)*16);
		memset(&this->countsAlleles[0],   0, sizeof(U64)*3);
	}
	~TomahawkImportRLEHelper(){}

	U64& operator[](const U32& p){ return(this->countsGenotypes[p]); }

	inline void reset(){
		this->countsGenotypes[0] = 0;
		this->countsGenotypes[1] = 0;
		this->countsGenotypes[4] = 0;
		this->countsGenotypes[5] = 0;
		this->countsAlleles[0]   = 0; // p
		this->countsAlleles[1]   = 0; // q
		this->countsAlleles[2]   = 0; // missing
	}
	inline const bool& hasMissing(void) const{ return(this->missingValues); }

	inline const double calculateAF(void) const{
		return(this->countsAlleles[0]/((double)this->countsAlleles[0]+this->countsAlleles[1]));
	}

	void calculateMGF(void){
		// Find the largest non-missing value
		U64 curMax = 0;
		U64* target; // To compare pointer address
		if(this->countsGenotypes[0] > curMax){
			curMax = this->countsGenotypes[0];
			target = &this->countsGenotypes[0];
		}
		if(this->countsGenotypes[1] > curMax){
			curMax = this->countsGenotypes[1];
			target = &this->countsGenotypes[1];
		}
		if(this->countsGenotypes[4] > curMax){
			curMax = this->countsGenotypes[4];
			target = &this->countsGenotypes[4];
		}
		if(this->countsGenotypes[5] > curMax){
			curMax = this->countsGenotypes[5];
			target = &this->countsGenotypes[5];
		}

		// Find next largest non-missing value
		U64 curMax2 = 0;
		if(this->countsGenotypes[0] >= curMax2 && &this->countsGenotypes[0] != target)
			curMax2 = this->countsGenotypes[0];
		if(this->countsGenotypes[1] >= curMax2 && &this->countsGenotypes[1] != target)
			curMax2 = this->countsGenotypes[1];
		if(this->countsGenotypes[4] >= curMax2 && &this->countsGenotypes[4] != target)
			curMax2 = this->countsGenotypes[4];
		if(this->countsGenotypes[5] >= curMax2 && &this->countsGenotypes[5] != target)
			curMax2 = this->countsGenotypes[5];

		this->MGF = (double)curMax2/this->expectedSamples;
	}

	bool determinePhase(const char& separator){
		if(separator == '/'){
			this->phased = false;
			return true;
		} else if(separator == '|'){
			this->phased = true;
			return true;
		} else return false;
	}

	/*
	// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
	// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
	// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
	//
	// Written by Jan Wigginton
	// Modified to use Tomahawk data
	*/
	void calculateHardyWeinberg(void){
		U64 obs_hets = this->countsGenotypes[1] + this->countsGenotypes[4];
		U64 obs_hom1 = this->countsGenotypes[0];
		U64 obs_hom2 = this->countsGenotypes[5];

		U64 obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
		U64 obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

		int64_t rare_copies = 2 * obs_homr + obs_hets;
		int64_t genotypes   = obs_hets + obs_homc + obs_homr;

		double* het_probs = new double[rare_copies + 1];

		int64_t i;
		for (i = 0; i <= rare_copies; ++i)
			het_probs[i] = 0.0;

		/* start at midpoint */
		int64_t mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

		/* check to ensure that midpoint and rare alleles have same parity */
		if ((rare_copies & 1) ^ (mid & 1))
			++mid;

		int64_t curr_hets = mid;
		int64_t curr_homr = (rare_copies - mid) / 2;
		int64_t curr_homc = genotypes - curr_hets - curr_homr;

		het_probs[mid] = 1.0;
		double sum = het_probs[mid];
		for (curr_hets = mid; curr_hets > 1; curr_hets -= 2){
			het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
							   / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
			sum += het_probs[curr_hets - 2];

			/* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
			++curr_homr;
			++curr_homc;
		}

		curr_hets = mid;
		curr_homr = (rare_copies - mid) / 2;
		curr_homc = genotypes - curr_hets - curr_homr;
		for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2){
			het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
							/((curr_hets + 2.0) * (curr_hets + 1.0));
			sum += het_probs[curr_hets + 2];

			/* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
			--curr_homr;
			--curr_homc;
		}

		for (i = 0; i <= rare_copies; i++)
			het_probs[i] /= sum;

		double p_hwe = 0.0;
		/*  p-value calculation for p_hwe  */
		for (i = 0; i <= rare_copies; i++){
			if (het_probs[i] > het_probs[obs_hets])
				continue;

			p_hwe += het_probs[i];
		}

		p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

		delete [] het_probs;

		this->HWE_P = p_hwe;
	}

	friend std::ostream& operator<<(std::ostream& os, const self_type& self){
		os << self.countsGenotypes[0] << '\t' << self.countsGenotypes[1] << '\t' << self.countsGenotypes[4] << '\t' << self.countsGenotypes[5] << '\t' << self.MGF;
		return(os);
	}

	inline const U64 countAlleles(void) const{ return(this->countsAlleles[0] + this->countsAlleles[1] + this->countsAlleles[2]); }

	U64   countsGenotypes[16];
	U64   countsAlleles[3];
	double MAF;
	double MGF;
	double HWE_P;
	bool  missingValues;
	bool  phased;
	const U64  expectedSamples;
};

class GenotypeEncoder {
	typedef GenotypeEncoder self_type;
	typedef MetaEntry       meta_entry_type;
	typedef ImporterFilters filter_type;
	typedef bool (self_type::*rleFunction)(const vcf::VCFLine& line, meta_entry_type& meta, io::BasicBuffer& runs); // Type cast pointer to function
	typedef bool (self_type::*bcfFunction)(const bcf::BCFEntry& line, meta_entry_type& meta, io::BasicBuffer& runs); // Type cast pointer to function

	typedef TomahawkImportRLEHelper helper_type;

public:
	GenotypeEncoder(const U64 samples, const filter_type& filters) :
		n_samples(samples),
		encode(nullptr),
		encodeComplex(nullptr),
		encodeBCF(nullptr),
		bit_width(0),
		shiftSize(0),
		helper(samples),
		savings(0),
		filters(filters)
	{
	}

	~GenotypeEncoder(){
	}

	void DetermineBitWidth(void){
		if(this->n_samples <= constants::UPPER_LIMIT_SAMPLES_8B - 1){
			if(!SILENT){
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << constants::UPPER_LIMIT_SAMPLES_8B << "..." << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Using 8-bit width..." << std::endl;
			}
			this->encode = &self_type::RunLengthEncodeSimple<BYTE>;
			this->encodeComplex = &self_type::RunLengthEncodeComplex<BYTE>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<BYTE>;
			this->shiftSize = sizeof(BYTE)*8 - constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width = sizeof(BYTE);
		} else if(this->n_samples <= constants::UPPER_LIMIT_SAMPLES_16B - 1){
			if(!SILENT){
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << constants::UPPER_LIMIT_SAMPLES_16B << "..." << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Using 16-bit width..." << std::endl;
			}
			this->encode = &self_type::RunLengthEncodeSimple<U16>;
			this->encodeComplex = &self_type::RunLengthEncodeComplex<U16>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<U16>;
			this->shiftSize = sizeof(U16)*8 - constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width = sizeof(U16);
		} else if(this->n_samples <= constants::UPPER_LIMIT_SAMPLES_32B - 1){
			if(!SILENT){
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << constants::UPPER_LIMIT_SAMPLES_32B << "..." << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Using 32-bit width..." << std::endl;
			}
			this->encode = &self_type::RunLengthEncodeSimple<U32>;
			this->encodeComplex = &self_type::RunLengthEncodeComplex<U32>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<U32>;
			this->shiftSize = sizeof(U32)*8 - constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width = sizeof(U32);
		} else {
			if(!SILENT){
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << constants::UPPER_LIMIT_SAMPLES_32B << "... Skip" << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << constants::UPPER_LIMIT_SAMPLES_64B << "..." << std::endl;
				std::cerr << helpers::timestamp("LOG", "RLE") << "Using 64-bit width..." << std::endl;
			}
			this->encode = &self_type::RunLengthEncodeSimple<U64>;
			this->encodeComplex = &self_type::RunLengthEncodeComplex<U64>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<U64>;
			this->shiftSize = sizeof(U64)*8 - constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width = sizeof(U64);
		}
	}

	inline bool RunLengthEncode(const vcf::VCFLine& line, meta_entry_type& meta, io::BasicBuffer& runs){
		if(!line.getComplex())
			return((*this.*encode)(line, meta, runs));
		else
			return((*this.*encodeComplex)(line, meta, runs));
	}

	inline bool RunLengthEncode(const bcf::BCFEntry& line, meta_entry_type& meta, io::BasicBuffer& runs){
		return((*this.*encodeBCF)(line, meta, runs));
	}

	inline const BYTE& getBitWidth(void) const{ return this->bit_width; }

private:
	template <class T> bool RunLengthEncodeSimple (const vcf::VCFLine& line, meta_entry_type& meta, io::BasicBuffer& runs);
	template <class T> bool RunLengthEncodeComplex(const vcf::VCFLine& line, meta_entry_type& meta, io::BasicBuffer& runs);
	template <class T> bool RunLengthEncodeBCF(const bcf::BCFEntry& line, meta_entry_type& meta, io::BasicBuffer& runs);

private:
	U64         n_samples;
	rleFunction encode;        // encoding function
	rleFunction encodeComplex; // encoding function
	bcfFunction encodeBCF;
	BYTE        bit_width;
	BYTE        shiftSize;     // bit shift size
	helper_type helper;
	const filter_type& filters;

public:
	U64 savings;
};

template <class T>
bool GenotypeEncoder::RunLengthEncodeBCF(const bcf::BCFEntry& line, meta_entry_type& meta, io::BasicBuffer& runs){
	meta.position = (U32)line.body->POS + 1;
	meta.ref_alt  = line.ref_alt;


	this->helper.reset();
	U32 internal_pos = line.formatID[0].l_offset;

	const std::vector<BYTE> swap = {0, 1, 2, 1, 0, 2};
	BYTE add = 0;

	// run pre-check
	{
		for(U32 i = 0; i < this->n_samples * 2; i += 2){
			const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
			const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
			BYTE packed_internal = (bcf::BCF_UNPACK_GENOTYPE(fmt_type_value2) << 2) | bcf::BCF_UNPACK_GENOTYPE(fmt_type_value1);
			this->helper[packed_internal] += 1;
			this->helper.countsAlleles[bcf::BCF_UNPACK_GENOTYPE(fmt_type_value2)] += 1;
			this->helper.countsAlleles[bcf::BCF_UNPACK_GENOTYPE(fmt_type_value1)] += 1;
		}

		// Univariate for reference allele
		if(this->filters.dropUnivariantRef && (this->helper.countsAlleles[0] == this->helper.countsAlleles[0] + this->helper.countsAlleles[1])){
			//std::cerr << "all reference: " << this->helper.countsAlleles[0] << "/" << this->helper.countsAlleles[1] << std::endl;
			return false;
		}

		// Univariate for alternative allele
		if(this->filters.dropUnivariantAlt && (this->helper.countsAlleles[1] == this->helper.countsAlleles[0] + this->helper.countsAlleles[1])){
			//std::cerr << "all alt: " << this->helper.countsAlleles[0] << "/" << this->helper.countsAlleles[1] << std::endl;
			return false;
		}

		// Flip reference and alternative allele such that 0 is the major allele
		if(this->filters.flipMajorMinor && this->helper.calculateAF() < 0.5){ // reference allele frequency
			//std::cerr << "alt > ref: " << this->helper.countsAlleles[0] << "/" << this->helper.countsAlleles[1] << std::endl;
			add = 3;
		}
		this->helper.reset();
	}
	internal_pos = line.formatID[0].l_offset;

	U64 sumLength = 0;
	T length = 1;
	T __dump = 0;
	T n_runs = 0;

	const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
	const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
	BYTE packed = (swap[bcf::BCF_UNPACK_GENOTYPE(fmt_type_value2) + add] << 2) | swap[bcf::BCF_UNPACK_GENOTYPE(fmt_type_value1) + add];

	this->helper.phased = fmt_type_value2 & 1; // MSB contains phasing information

	for(U32 i = 2; i < this->n_samples * 2; i += 2){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
		BYTE packed_internal = (swap[bcf::BCF_UNPACK_GENOTYPE(fmt_type_value2) + add] << 2) | swap[bcf::BCF_UNPACK_GENOTYPE(fmt_type_value1) + add];

		if(packed != packed_internal){
			__dump =  (length & (((T)1 << this->shiftSize) - 1)) << constants::TOMAHAWK_SNP_PACK_WIDTH;
			__dump ^= (packed & ((1 << constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
			runs += __dump;

			this->helper[packed] += length;
			this->helper.countsAlleles[packed >> 2] += length;
			this->helper.countsAlleles[packed & 3]  += length;

			sumLength += length;
			length = 1;
			packed = packed_internal;
			++n_runs;
			continue;
		}
		++length;
	}
	__dump =  (length & (((T)1 << this->shiftSize) - 1)) << constants::TOMAHAWK_SNP_PACK_WIDTH;
	__dump ^= (packed & ((1 << constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
	runs += __dump;
	++n_runs;

	this->helper[packed] += length;
	this->helper.countsAlleles[packed >> 2] += length;
	this->helper.countsAlleles[packed & 3]  += length;

	sumLength += length;
	assert(sumLength == this->n_samples);

	this->helper.calculateMGF();
	this->helper.calculateHardyWeinberg();

	meta.all_phased  = this->helper.phased;
	meta.has_missing = this->helper.missingValues;
	meta.AF          = this->helper.calculateAF();
	meta.HWE_P       = this->helper.HWE_P;
	meta.runs        = n_runs;

	this->helper.reset();
	return true;
}

template <class T>
bool GenotypeEncoder::RunLengthEncodeSimple(const vcf::VCFLine& line, meta_entry_type& meta, io::BasicBuffer& runs){
	return false;
}


template <class T>
bool GenotypeEncoder::RunLengthEncodeComplex(const vcf::VCFLine& line, meta_entry_type& meta, io::BasicBuffer& runs){
	return false;
}

}
}

#endif /* TOMAHAWKIMPORTRLE_H_ */
