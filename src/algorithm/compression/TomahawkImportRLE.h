#ifndef TOMAHAWKIMPORTRLE_H_
#define TOMAHAWKIMPORTRLE_H_

#include <algorithm>
#include <bitset>

#include "../../math/FisherMath.h"
#include "../../tomahawk/base/TomahawkEntryMeta.h"
#include "../../io/bcf/BCFReader.h"
#include "RunLengthEncoding.h"

namespace Tomahawk{
namespace Algorithm{

// Todo: These should all be moved to VCF or BCF entry class
struct TomahawkImportRLEHelper{
	typedef TomahawkImportRLEHelper self_type;

	TomahawkImportRLEHelper(const U64 expectedSamples) :
		MAF(0),
		MGF(0),
		HWE_P(0),
		missingValues(0),
		phased(false),
		expectedSamples(expectedSamples),
		fisherTable(1)
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

	double calculateMAF(void){
		if(this->countsAlleles[0] > this->countsAlleles[1])
			return(this->countsAlleles[1]/((double)this->countsAlleles[0]+this->countsAlleles[1]));
		else
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

	U64 countsGenotypes[16];
	U64 countsAlleles[3];
	float MAF;
	float MGF;
	float HWE_P;
	bool missingValues;
	bool phased;
	const U64 expectedSamples;
	FisherMath fisherTable;
};

class TomahawkImportRLE {
	typedef TomahawkImportRLE self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef VCF::VCFLine vcf_type;
	typedef BCF::BCFEntry bcf_type;
	typedef Support::TomahawkEntryMetaBase meta_base_type;
	typedef bool (Tomahawk::Algorithm::TomahawkImportRLE::*rleFunction)(const vcf_type& line, buffer_type& meta, buffer_type& runs); // Type cast pointer to function
	typedef bool (Tomahawk::Algorithm::TomahawkImportRLE::*bcfFunction)(const bcf_type& line, meta_base_type& meta_base, buffer_type& meta, buffer_type& runs); // Type cast pointer to function
	typedef TomahawkImportRLEHelper helper_type;

public:
	TomahawkImportRLE(const U64 samples) :
		n_samples(samples),
		encode(nullptr),
		encodeComplex(nullptr),
		encodeBCF(nullptr),
		bit_width(0),
		shiftSize(0),
		helper(samples),
		savings(0)
	{
	}

	~TomahawkImportRLE(){
	}

	void DetermineBitWidth(void){
		if(this->n_samples <= Constants::UPPER_LIMIT_SAMPLES_8B){
			if(!SILENT){
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << Constants::UPPER_LIMIT_SAMPLES_8B << "..." << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 8-bit width..." << std::endl;
			}
			this->encode = &self_type::RunLengthEncodeSimple<BYTE>;
			this->encodeComplex = &self_type::RunLengthEncodeComplex<BYTE>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<BYTE>;
			this->shiftSize = sizeof(BYTE)*8 - Constants::TOMAHAWK_SHIFT_SIZE;
			this->bit_width = sizeof(BYTE);
		} else if(this->n_samples <= Constants::UPPER_LIMIT_SAMPLES_16B){
			if(!SILENT){
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << Constants::UPPER_LIMIT_SAMPLES_16B << "..." << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 16-bit width..." << std::endl;
			}
			this->encode = &self_type::RunLengthEncodeSimple<U16>;
			this->encodeComplex = &self_type::RunLengthEncodeComplex<U16>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<U16>;
			this->shiftSize = sizeof(U16)*8 - Constants::TOMAHAWK_SHIFT_SIZE;
			this->bit_width = sizeof(U16);
		} else if(this->n_samples <= Constants::UPPER_LIMIT_SAMPLES_32B){
			if(!SILENT){
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << Constants::UPPER_LIMIT_SAMPLES_32B << "..." << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 32-bit width..." << std::endl;
			}
			this->encode = &self_type::RunLengthEncodeSimple<U32>;
			this->encodeComplex = &self_type::RunLengthEncodeComplex<U32>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<U32>;
			this->shiftSize = sizeof(U32)*8 - Constants::TOMAHAWK_SHIFT_SIZE;
			this->bit_width = sizeof(U32);
		} else {
			if(!SILENT){
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_32B << "... Skip" << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << Constants::UPPER_LIMIT_SAMPLES_64B << "..." << std::endl;
				std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 64-bit width..." << std::endl;
			}
			this->encode = &self_type::RunLengthEncodeSimple<U64>;
			this->encodeComplex = &self_type::RunLengthEncodeComplex<U64>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<U64>;
			this->shiftSize = sizeof(U64)*8 - Constants::TOMAHAWK_SHIFT_SIZE;
			this->bit_width = sizeof(U64);
		}
	}

	inline bool RunLengthEncode(const vcf_type& line, buffer_type& meta, buffer_type& runs){
		if(!line.getComplex())
			return((*this.*encode)(line, meta, runs));
		else
			return((*this.*encodeComplex)(line, meta, runs));
	}

	inline bool RunLengthEncode(const bcf_type& line, meta_base_type& meta_base, buffer_type& meta, buffer_type& runs){
		return((*this.*encodeBCF)(line, meta_base, meta, runs));
	}

	inline const BYTE& getBitWidth(void) const{ return this->bit_width; }

private:
	template <class T> bool RunLengthEncodeSimple (const vcf_type& line, buffer_type& meta, buffer_type& runs);
	template <class T> bool RunLengthEncodeComplex(const vcf_type& line, buffer_type& meta, buffer_type& runs);
	template <class T> bool RunLengthEncodeBCF(const bcf_type& line, meta_base_type& meta_base, buffer_type& meta, buffer_type& runs);

private:
	U64 n_samples;
	rleFunction encode;			// encoding function
	rleFunction encodeComplex;	// encoding function
	bcfFunction encodeBCF;
	BYTE bit_width;
	BYTE shiftSize;				// bit shift size
	helper_type helper;

public:
	U64 savings;
};

template <class T>
bool TomahawkImportRLE::RunLengthEncodeBCF(const bcf_type& line, meta_base_type& meta_base, buffer_type& meta, buffer_type& runs){
	// Update basic values for a meta entry
	meta_base.position = (U32)line.body->POS + 1; // Base-1
	meta_base.ref_alt = line.ref_alt;
	meta_base.biallelic = true;
	meta_base.simple = line.isSimple();
	T n_runs = 0;

	// If the number of alleles == 2
	if(line.body->n_allele == 2){
		U32 internal_pos = line.p_genotypes; // virtual byte offset of genotype start
		U64 sumLength = 0;
		T length = 1;
		T RLE = 0;

		const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
		BYTE packed = (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) << 3) | (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) << 1) | (fmt_type_value2 & 1);

		// MSB contains phasing information
		this->helper.phased = (packed & 1);

		for(U32 i = 2; i < this->n_samples * 2; i += 2){
			const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
			const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
			const BYTE packed_internal = (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) << 3) | (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) << 1) | (fmt_type_value2 & 1);

			if(packed != packed_internal){
				// Prepare RLE
				RLE = length;
				RLE <<= Constants::TOMAHAWK_SHIFT_SIZE;
				RLE |= packed;

				// Set meta phased flag bit
				if((packed & 1) != 1) this->helper.phased = false;

				// Push RLE to buffer
				runs += RLE;

				// Update counts
				this->helper[packed >> 1] += length;
				this->helper.countsAlleles[packed >> 3] += length;
				this->helper.countsAlleles[(packed >> 1) & 3]  += length;

				// Reset and update
				sumLength += length;
				length = 1;
				packed = packed_internal;
				++n_runs;
				continue;
			}
			++length;
		}
		// Last entry
		// Prepare RLE
		RLE = length;
		RLE <<= Constants::TOMAHAWK_SHIFT_SIZE;
		RLE |= packed;

		// Set meta phased flag bit
		if((packed & 1) != 1) this->helper.phased = false;

		// Push RLE to buffer
		runs += RLE;
		++n_runs;

		// Update counts
		this->helper[packed >> 1] += length;
		this->helper.countsAlleles[packed >> 3] += length;
		this->helper.countsAlleles[(packed >> 1) & 3]  += length;

		// Reset and update
		sumLength += length;
		assert(sumLength == this->n_samples);

		// Calculate basic stats
		this->helper.calculateMGF();
		this->helper.calculateHardyWeinberg();
	}
	// If number of alleles != 2
	// Revert back to no-encoding
	else {
		// We use ceil(log2(n_alleles + 1)) bits for each allele
		// the + 1 represent the missing case
		// MSB first bit encode for phasing between the diploid alleles
		// The word size is then (2*bits + 1) / 8 where + 1 represent the phasing
		// information
		const U32 shift_size = ceil(log2(double(line.body->n_allele) + 1));
		const BYTE bit_width = 2*shift_size + 1;
		const BYTE word_size = ceil(float(bit_width) / 8);
		meta_base.biallelic = false;

		// Virtual byte offset into start of genotypes
		// in BCF entry
		U32 internal_pos = line.p_genotypes;
		n_runs = this->n_samples;

		// Todo: move to templated function
		// If word size is BYTE
		if(word_size == 1){
			for(U32 i = 0; i < this->n_samples * 2; i += 2){
				const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
				const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
				const BYTE packed = (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) << (shift_size + 1)) | (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) << 1) | (fmt_type_value2 & 1);
				runs += packed;
			}
		}
		// If word_size is U16
		else {
			std::cerr << Helpers::timestamp("ERROR", "RLE") << "Not implemented word-size > 1..." << std::endl;
			exit(1);
		}

		this->helper.HWE_P = 1;
		this->helper.MGF = 0;
	}

	// See meta for more information
	meta_base.phased = this->helper.phased;
	meta_base.missing = this->helper.missingValues;
	meta_base.MGF = this->helper.MGF;
	meta_base.HWE_P = this->helper.HWE_P;
	meta_base.virtual_offset = runs.pointer; // position at end
	meta += meta_base;
	meta += n_runs;

	// Reset and recycle helper
	this->helper.reset();
	return true;
}

template <class T>
bool TomahawkImportRLE::RunLengthEncodeSimple(const vcf_type& line, buffer_type& meta, buffer_type& runs){
	meta += line.position;
	meta += line.ref_alt;

	///////////////////////////////
	// Encoding:
	// First 8|T| - TOMAHAWK_SNP_PACK_WIDTH bits encode the run length
	// remaining TOMAHAWK_SNP_PACK_WIDTH bits encode
	// TOMAHAWK_ALLELE_PACK_WIDTH bits of snpA and TOMAHAWK_ALLELE_PACK_WIDTH bits of snpB
	///////////////////////////////
	T run_length = 1;

	// ASCII value for '.' is 46
	// Therefore:
	// . - 46 = 0
	// 0 - 46 = 2
	// 1 - 46 = 3
	//
	// Remap:
	// 0 -> 2 (missing)
	// 2 -> 0 (1)
	// 3 -> 1 (0)
	//
	// Genotypes are thus:
	// 0/0 -> 0000b = 0
	// 0/1 -> 0001b = 1
	// 1/0 -> 0100b = 4
	// 1/1 -> 0101b = 5
	// ...
	// mixed values for missing
	// ....
	// largest value:
	// ./. -> 0101b = 10

	// Determine phase of the variant
	// Sets the phase of all genotypes in this variant to whatever the first one is
	this->helper.determinePhase(line.simple_[0].separator);

	BYTE type =  Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[0].snpA - 46] << Constants::TOMAHAWK_ALLELE_PACK_WIDTH;
		 type ^= Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[0].snpB - 46] << 0;
	BYTE curType = type;
	T __dump = 0;
	T total_samples = 0;
	T runsCount = 0;
	//U32 startOffset = runs.pointer;

	// Encode
	for(U32 i = 1; i < this->n_samples; ++i){
		curType =  Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[i].snpA - 46] << Constants::TOMAHAWK_ALLELE_PACK_WIDTH;
		curType ^= Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[i].snpB - 46] << 0;

		// Current type is different from previous type
		if(curType != type){
			__dump =  (run_length & (((T)1 << this->shiftSize) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
			__dump ^= (type & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));


			total_samples += run_length;
			//runs.pointer += PACK3(curType, &runs.data[runs.pointer], run_length);
			runs += __dump;

			this->helper[type] += run_length;
			this->helper.countsAlleles[type >> 2] += run_length;
			this->helper.countsAlleles[type & 3]  += run_length;

			run_length = 1;
			type = curType;
			++runsCount;

		} else ++run_length;
	}

	// Encode final
	__dump =  (run_length & (((T)1 << this->shiftSize) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
	__dump ^= (type & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
	runs += __dump;
	total_samples += run_length;
	//runs.pointer += PACK3(curType, &runs.data[runs.pointer], run_length);
	this->helper[type] += run_length;
	this->helper.countsAlleles[type >> 2] += run_length;
	this->helper.countsAlleles[type & 3]  += run_length;
	++runsCount;

	if(total_samples != this->n_samples){
		std::cerr << Helpers::timestamp("ERROR", "RLE") << "Sum of run lengths does not equal number of samples: " << total_samples << "/" << this->n_samples << std::endl;
		exit(1);
	}

	this->helper.calculateMGF();
	this->helper.calculateHardyWeinberg();

	// Position
	U32& position = *reinterpret_cast<U32*>(&meta[meta.pointer - 5]);
	position <<= 2;
	position |= this->helper.phased << 1;
	position |= this->helper.missingValues << 0;
	meta += this->helper.MGF;
	meta += this->helper.HWE_P;
	meta += runsCount;

	this->helper.reset();
	return true;
}


template <class T>
bool TomahawkImportRLE::RunLengthEncodeComplex(const vcf_type& line, buffer_type& meta, buffer_type& runs){
	meta += line.position;
	meta += line.ref_alt;

	///////////////////////////////
	// Encoding:
	// First 8|T| - TOMAHAWK_SNP_PACK_WIDTH bits encode the run length
	// remaining TOMAHAWK_SNP_PACK_WIDTH bits encode
	// TOMAHAWK_ALLELE_PACK_WIDTH bits of snpA and TOMAHAWK_ALLELE_PACK_WIDTH bits of snpB
	///////////////////////////////
	T run_length = 1;

	// ASCII value for '.' is 46
	// Therefore:
	// . - 46 = 0
	// 0 - 46 = 2
	// 1 - 46 = 3
	//
	// Remap:
	// 0 -> 2 (missing)
	// 2 -> 0 (1)
	// 3 -> 1 (0)
	//
	// Genotypes are thus:
	// 0/0 -> 0000b = 0
	// 0/1 -> 0001b = 1
	// 1/0 -> 0100b = 4
	// 1/1 -> 0101b = 5
	// ...
	// mixed values for missing
	// ....
	// largest value:
	// ./. -> 0101b = 10

	// Determine phase of the variant
	// Sets the phase of all genotypes in this variant to whatever the first one is
	this->helper.determinePhase(line.complex_[0]->separator);

	BYTE type =  Constants::TOMAHAWK_ALLELE_LOOKUP[line.complex_[0]->snpA - 46] << Constants::TOMAHAWK_ALLELE_PACK_WIDTH;
		 type ^= Constants::TOMAHAWK_ALLELE_LOOKUP[line.complex_[0]->snpB - 46] << 0;
	BYTE curType = type;
	T __dump = 0;
	T total_samples = 0;
	T runsCount = 0;

	// Encode
	for(U32 i = 1; i < this->n_samples; ++i){
		curType =  Constants::TOMAHAWK_ALLELE_LOOKUP[line.complex_[i]->snpA - 46] << Constants::TOMAHAWK_ALLELE_PACK_WIDTH;
		curType ^= Constants::TOMAHAWK_ALLELE_LOOKUP[line.complex_[i]->snpB - 46] << 0;

		// Current type is different from previous type
		if(curType != type){
			__dump =  (run_length & (((T)1 << this->shiftSize) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
			__dump ^= (type & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
			runs += __dump;
			this->helper[type] += run_length;
			this->helper.countsAlleles[type >> 2] += run_length;
			this->helper.countsAlleles[type & 3]  += run_length;

			//std::cerr << run_length << '|' << (int)type << '\t';

			total_samples += run_length;
			run_length = 1;
			type = curType;
			++runsCount;

		} else ++run_length;
	}

	// Encode final
	__dump =  (run_length & (((T)1 << this->shiftSize) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
	__dump ^= (curType & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
	runs += __dump;
	this->helper[type] += run_length;
	this->helper.countsAlleles[type >> 2] += run_length;
	this->helper.countsAlleles[type & 3]  += run_length;
	//std::cerr << run_length << '|' << (int)type << std::endl;
	++runsCount;

	total_samples += run_length;

	if(total_samples != this->n_samples){
		std::cerr << Helpers::timestamp("ERROR", "RLE") << "Sum of run lengths does not equal number of samples: " << total_samples << "/" << this->n_samples << std::endl;
		exit(1);
	}

	this->helper.calculateMGF();
	this->helper.calculateHardyWeinberg();

	// Position
	U32& position = *reinterpret_cast<U32*>(&meta[meta.pointer - 5]);
	position <<= 2;
	position |= this->helper.phased << 1;
	position |= this->helper.missingValues << 0;
	meta += this->helper.MGF;
	meta += this->helper.HWE_P;
	meta += runsCount;

	this->helper.reset();
	return true;
}

}
}

#endif /* TOMAHAWKIMPORTRLE_H_ */
