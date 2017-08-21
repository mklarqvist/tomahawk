#ifndef TOMAHAWKIMPORTRLE_H_
#define TOMAHAWKIMPORTRLE_H_

#include <algorithm>
#include <bitset>

#include "../../math/FisherTest.h"
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

	void determinePhase(const char& separator){
		switch(separator){
		case '/': this->phased = false; break;
		case '|': this->phased = true; break;
		default:
			std::cerr << "ERROR" << std::endl;
			exit(1);
			break;
		}
	}

	void calculateHardyWeinberg(void){
		// Total number of non-missing genotypes
		const U64 totalValidGenotypes = this->countsGenotypes[0] + this->countsGenotypes[1] + this->countsGenotypes[4] + this->countsGenotypes[5];

		// If total valid genotypes is not equal to all the individuals in a line
		// trigger missing values flag
		if(totalValidGenotypes != this->expectedSamples)
			this->missingValues = true;

		// Calculate P and Q frequencies from genotypes
		//
		// 0 -> Number of ref-ref
		// 1 -> Number of alt-ref
		// 4 -> Number of ref-alt
		// 5 -> Number of alt-alt
		//
		// P is therefore: 2*{0} + {1} + {4} / 2*{0,1,4,5}
		// Q is therefore: 2*{5} + {1} + {4} / 2*{0,1,4,5}
		const double p = ((double)2*this->countsGenotypes[0]+this->countsGenotypes[1]+this->countsGenotypes[4])/(2*totalValidGenotypes);
		const double q = ((double)this->countsGenotypes[1]+this->countsGenotypes[4]+2*this->countsGenotypes[5])/(2*totalValidGenotypes);
		const double pp = p*p*totalValidGenotypes; // p^2
		const double pq = 2*p*q*totalValidGenotypes; // 2pq
		const double qq = q*q*totalValidGenotypes; // q^2
		double ppCV = pow(this->countsGenotypes[0] - pp,2)/pp;
		double pqCV = pow((this->countsGenotypes[1] + this->countsGenotypes[4]) - pq,2)/pq;
		double qqCV = pow(this->countsGenotypes[5] - qq,2)/qq;
		if(pp < 1) ppCV = 0;
		if(pq < 1) pqCV = 0;
		if(qq < 1) qqCV = 0;

		const double CV = ppCV + pqCV + qqCV;
		this->HWE_P = this->fisherTable.chisqr(1, CV);
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
	FisherTest fisherTable;
};

class TomahawkImportRLE{
	typedef TomahawkImportRLE self_type;
	typedef void (Tomahawk::Algorithm::TomahawkImportRLE::*rleFunction)(const VCF::VCFLine& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs); // Type cast pointer to function
	typedef void (Tomahawk::Algorithm::TomahawkImportRLE::*bcfFunction)(const BCF::BCFEntry& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs); // Type cast pointer to function

	typedef TomahawkImportRLEHelper helper_type;

public:
	TomahawkImportRLE(const U64 samples) :
		n_samples(samples),
		encode_(nullptr),
		encodeComplex_(nullptr),
		encodeBCF(nullptr),
		bit_width_(0),
		shiftSize_(0),
		helper_(samples),
		savings(0)
	{
	}

	~TomahawkImportRLE(){
	}

	void DetermineBitWidth(void){
		if(this->n_samples <= Constants::UPPER_LIMIT_SAMPLES_8B - 1){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << Constants::UPPER_LIMIT_SAMPLES_8B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 8-bit width..." << std::endl;
			this->encode_ = &self_type::RunLengthEncodeSimple<BYTE>;
			this->encodeComplex_ = &self_type::RunLengthEncodeComplex<BYTE>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<BYTE>;
			this->shiftSize_ = sizeof(BYTE)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width_ = sizeof(BYTE);
		} else if(this->n_samples <= Constants::UPPER_LIMIT_SAMPLES_16B - 1){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << Constants::UPPER_LIMIT_SAMPLES_16B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 16-bit width..." << std::endl;
			this->encode_ = &self_type::RunLengthEncodeSimple<U16>;
			this->encodeComplex_ = &self_type::RunLengthEncodeComplex<U16>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<U16>;
			this->shiftSize_ = sizeof(U16)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width_ = sizeof(U16);
		} else if(this->n_samples <= Constants::UPPER_LIMIT_SAMPLES_32B - 1){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << Constants::UPPER_LIMIT_SAMPLES_32B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 32-bit width..." << std::endl;
			this->encode_ = &self_type::RunLengthEncodeSimple<U32>;
			this->encodeComplex_ = &self_type::RunLengthEncodeComplex<U32>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<U32>;
			this->shiftSize_ = sizeof(U32)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width_ = sizeof(U32);
		} else {
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " > " << Constants::UPPER_LIMIT_SAMPLES_32B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->n_samples << " < " << Constants::UPPER_LIMIT_SAMPLES_64B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 64-bit width..." << std::endl;
			this->encode_ = &self_type::RunLengthEncodeSimple<U64>;
			this->encodeComplex_ = &self_type::RunLengthEncodeComplex<U64>;
			this->encodeBCF = &self_type::RunLengthEncodeBCF<U64>;
			this->shiftSize_ = sizeof(U64)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width_ = sizeof(U64);
		}
	}

	inline void RunLengthEncode(const VCF::VCFLine& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs){
		if(!line.getComplex())
			(*this.*encode_)(line, meta, runs);
		else
			(*this.*encodeComplex_)(line, meta, runs);
	}

	inline void RunLengthEncode(const BCF::BCFEntry& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs){
		return((*this.*encodeBCF)(line, meta, runs));
	}

	template <class T>
	void RunLengthEncodeBCF(const BCF::BCFEntry& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs){
		U32 internal_pos = line.p_genotypes;
		T length = 1;
		U64 sumLength = 0;
		T __dump = 0;
		T n_runs = 0;

		const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
		BYTE packed = (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) << 2) | BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2);

		for(U32 i = 2; i < this->n_samples * 2; i+=2){
			const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
			const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&line.data[internal_pos++]);
			const BYTE packed_internal = (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) << 2) | BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2);
			if(packed != packed_internal){
				__dump =  (length & (((T)1 << this->shiftSize_) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
				__dump ^= (packed & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
				runs += __dump;

				//std::cerr << length << '|' << (int)packed << std::endl;

				this->helper_[packed] += length;
				this->helper_.countsAlleles[packed >> 2] += length;
				this->helper_.countsAlleles[packed & 3]  += length;

				sumLength += length;
				length = 1;
				packed = packed_internal;
				++n_runs;
				continue;
			}
			++length;
		}
		__dump =  (length & (((T)1 << this->shiftSize_) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
		__dump ^= (packed & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
		runs += __dump;
		++n_runs;

		this->helper_[packed] += length;
		this->helper_.countsAlleles[packed >> 2] += length;
		this->helper_.countsAlleles[packed & 3]  += length;

		sumLength += length;
		assert(sumLength == this->n_samples);
		assert(internal_pos == line.size());

		this->helper_.calculateMGF();
		this->helper_.calculateHardyWeinberg();

		// Position
		U32& position = *reinterpret_cast<U32*>(&meta[meta.pointer - 5]);
		position <<= 2;
		position |= this->helper_.phased << 1;
		position |= this->helper_.missingValues << 0;
		meta += this->helper_.MGF;
		meta += this->helper_.HWE_P;
		meta += n_runs;

		this->helper_.reset();
	}

	inline const BYTE& getBitWidth(void) const{ return this->bit_width_; }

private:
	template <class T> void RunLengthEncodeSimple (const VCF::VCFLine& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs);
	template <class T> void RunLengthEncodeComplex(const VCF::VCFLine& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs);

private:
	U64 n_samples;
	rleFunction encode_;			// encoding function
	rleFunction encodeComplex_;		// encoding function
	bcfFunction encodeBCF;
	BYTE bit_width_;
	BYTE shiftSize_;				// bit shift size
	helper_type helper_;
public:
	U64 savings;
};

template <class T>
void TomahawkImportRLE::RunLengthEncodeSimple(const VCF::VCFLine& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs){
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
	this->helper_.determinePhase(line.simple_[0].separator);

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
			__dump =  (run_length & (((T)1 << this->shiftSize_) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
			__dump ^= (type & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));


			total_samples += run_length;
			//runs.pointer += PACK3(curType, &runs.data[runs.pointer], run_length);
			runs += __dump;

			this->helper_[type] += run_length;
			this->helper_.countsAlleles[type >> 2] += run_length;
			this->helper_.countsAlleles[type & 3]  += run_length;


			run_length = 1;
			type = curType;
			++runsCount;

		} else ++run_length;
	}

	// Encode final
	__dump =  (run_length & (((T)1 << this->shiftSize_) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
	__dump ^= (type & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
	runs += __dump;
	total_samples += run_length;
	//runs.pointer += PACK3(curType, &runs.data[runs.pointer], run_length);
	this->helper_[type] += run_length;
	this->helper_.countsAlleles[type >> 2] += run_length;
	this->helper_.countsAlleles[type & 3]  += run_length;
	++runsCount;

	if(total_samples != this->n_samples){
		std::cerr << Helpers::timestamp("ERROR", "RLE") << "Sum of run lengths does not equal number of samples: " << total_samples << "/" << this->n_samples << std::endl;
		exit(1);
	}


	this->helper_.calculateMGF();
	this->helper_.calculateHardyWeinberg();

	//std::cerr << this->helper_.countsAlleles[0] << '\t' << this->helper_.countsAlleles[1] << '\t' << this->helper_.countsAlleles[2] << '\t' << this->helper_.countAlleles() << '\t' << this->helper_.calculateMAF() << '\t' << this->helper_.MGF << std::endl;


	// Position
	U32& position = *reinterpret_cast<U32*>(&meta[meta.pointer - 5]);
	position <<= 2;
	position |= this->helper_.phased << 1;
	position |= this->helper_.missingValues << 0;
	meta += this->helper_.MGF;
	meta += this->helper_.HWE_P;
	meta += runsCount;

	this->helper_.reset();
}


template <class T>
void TomahawkImportRLE::RunLengthEncodeComplex(const VCF::VCFLine& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs){
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
	this->helper_.determinePhase(line.complex_[0]->separator);

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
			__dump =  (run_length & (((T)1 << this->shiftSize_) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
			__dump ^= (type & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
			runs += __dump;
			this->helper_[type] += run_length;

			//std::cerr << run_length << '|' << (int)type << '\t';

			total_samples += run_length;
			run_length = 1;
			type = curType;
			++runsCount;

		} else ++run_length;
	}

	// Encode final
	__dump =  (run_length & (((T)1 << this->shiftSize_) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
	__dump ^= (curType & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
	runs += __dump;
	this->helper_[type] += run_length;
	//std::cerr << run_length << '|' << (int)type << std::endl;
	++runsCount;

	total_samples += run_length;

	if(total_samples != this->n_samples){
		std::cerr << Helpers::timestamp("ERROR", "RLE") << "Sum of run lengths does not equal number of samples: " << total_samples << "/" << this->n_samples << std::endl;
		exit(1);
	}

	this->helper_.calculateMGF();
	this->helper_.calculateHardyWeinberg();

	// Position
	U32& position = *reinterpret_cast<U32*>(&meta[meta.pointer - 5]);
	position <<= 2;
	position |= this->helper_.phased << 1;
	position |= this->helper_.missingValues << 0;

	meta += this->helper_.MGF;
	meta += this->helper_.HWE_P;
	meta += runsCount;

	this->helper_.reset();
}

}
}

#endif /* TOMAHAWKIMPORTRLE_H_ */
