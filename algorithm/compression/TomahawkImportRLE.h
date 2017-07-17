#ifndef TOMAHAWKIMPORTRLE_H_
#define TOMAHAWKIMPORTRLE_H_

#include <algorithm>
#include <bitset>

#include "../FisherTest.h"

namespace Tomahawk{
namespace Algorithm{


template <class T>
int static inline PACK3(const BYTE& ref, char* target, T& length){
	if(length <= 7){
		*target++ = 0 | ((ref & 15) << 3) | (length & 7); // highest bit is 0
		return 1;
	}

	char* target0 = target;
	*target++ = 128 | ((ref & 15) << 3) | (length & 7);
	length >>= 3;

	while(true){
		if(length <= 7){
			*target++ = 0 | (length & 127); // highest bit is 0
			length >>= 7;
			break;
		}

		*target++ = 128 | (length & 127);
		length >>= 7;
	}
	return(target - target0);
}

template <class T>
int static inline UNPACK3(char* target, T& length, BYTE& ref){
	if((*target & 128) == 0){
		length = *target & 7;
		ref = *target >> 3;
		return 1;
	}

	char* target0 = target;
	length = *target & 7;
	ref = (*target >> 3) & 15;
	++target;
	U32 offset = 3;

	while(true){
		length |= (*target & 127) << offset;
		offset += 7;

		if((*target & 128) == 0) break;
		++target;
	}

	return(target-target0);
}

struct TomahawkImportRLEHelper{
	typedef TomahawkImportRLEHelper self_type;

	TomahawkImportRLEHelper(const U64 expectedSamples) :
		MAF(0),
		HWE_P(0),
		missingValues(0),
		phased(false),
		expectedSamples(expectedSamples),
		fisherTable(1)
	{
		memset(&this->counts[0], 0, sizeof(U64)*16);
	}
	~TomahawkImportRLEHelper(){}

	U64& operator[](const U32& p){ return(this->counts[p]); }

	inline void reset(){ counts[0] = 0; counts[1] = 0; counts[4] = 0; counts[5] = 0; }
	inline const bool& hasMissing(void) const{ return(this->missingValues); }

	void calculateMAF(void){
		// Find the largest non-missing value
		U64 curMax = 0;
		U64* target; // To compare pointer address
		if(this->counts[0] > curMax){
			curMax = this->counts[0];
			target = &this->counts[0];
		}
		if(this->counts[1] > curMax){
			curMax = this->counts[1];
			target = &this->counts[1];
		}
		if(this->counts[4] > curMax){
			curMax = this->counts[4];
			target = &this->counts[4];
		}
		if(this->counts[5] > curMax){
			curMax = this->counts[5];
			target = &this->counts[5];
		}

		// Find next largest non-missing value
		U64 curMax2 = 0;
		if(this->counts[0] >= curMax2 && &this->counts[0] != target)
			curMax2 = this->counts[0];
		if(this->counts[1] >= curMax2 && &this->counts[1] != target)
			curMax2 = this->counts[1];
		if(this->counts[4] >= curMax2 && &this->counts[4] != target)
			curMax2 = this->counts[4];
		if(this->counts[5] >= curMax2 && &this->counts[5] != target)
			curMax2 = this->counts[5];

		this->MAF = (double)curMax2/this->expectedSamples;
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
		const U64 totalValidGenotypes = this->counts[0] + this->counts[1] + this->counts[4] + this->counts[5];

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
		const double p = ((double)2*this->counts[0]+this->counts[1]+this->counts[4])/(2*totalValidGenotypes);
		const double q = ((double)this->counts[1]+this->counts[4]+2*this->counts[5])/(2*totalValidGenotypes);
		const double pp = p*p*totalValidGenotypes;
		const double pq = 2*p*q*totalValidGenotypes;
		const double qq = q*q*totalValidGenotypes;
		double ppCV = pow(this->counts[0] - pp,2)/pp;
		double pqCV = pow((this->counts[1] + this->counts[4]) - pq,2)/pq;
		double qqCV = pow(this->counts[5] - qq,2)/qq;
		if(pp < 1) ppCV = 0;
		if(pq < 1) pqCV = 0;
		if(qq < 1) qqCV = 0;

		const double CV = ppCV + pqCV + qqCV;
		this->HWE_P = this->fisherTable.chisqr(1, CV);
	}

	friend std::ostream& operator<<(std::ostream& os, const self_type& self){
		os << self.counts[0] << '\t' << self.counts[1] << '\t' << self.counts[4] << '\t' << self.counts[5] << '\t' << self.MAF;
		return(os);
	}

	U64 counts[16];
	double MAF;
	double HWE_P;
	bool missingValues;
	bool phased;
	const U64 expectedSamples;
	FisherTest fisherTable;
};

class TomahawkImportRLE{
	typedef TomahawkImportRLE self_type;
	typedef void (Tomahawk::Algorithm::TomahawkImportRLE::*rleFunction)(const VCF::VCFLine& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs); // Type cast pointer to function
	typedef TomahawkImportRLEHelper helper_type;

public:
	TomahawkImportRLE(VCF::VCFHeader& header) :
		VCFheader_(header),
		encode_(nullptr),
		encodeComplex_(nullptr),
		bit_width_(0),
		shiftSize_(0),
		helper_(header.samples_),
		savings(0)
	{
	}

	~TomahawkImportRLE(){
	}

	void DetermineBitWidth(void){
		if(this->VCFheader_.size() <= Constants::UPPER_LIMIT_SAMPLES_8B - 1){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " < " << Constants::UPPER_LIMIT_SAMPLES_8B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 8-bit width..." << std::endl;
			this->encode_ = &TomahawkImportRLE::RunLengthEncodeSimple<BYTE>;
			this->encodeComplex_ = &TomahawkImportRLE::RunLengthEncodeComplex<BYTE>;
			this->shiftSize_ = sizeof(BYTE)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width_ = sizeof(BYTE);
		} else if(this->VCFheader_.size() <= Constants::UPPER_LIMIT_SAMPLES_16B - 1){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " < " << Constants::UPPER_LIMIT_SAMPLES_16B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 16-bit width..." << std::endl;
			this->encode_ = &TomahawkImportRLE::RunLengthEncodeSimple<U16>;
			this->encodeComplex_ = &TomahawkImportRLE::RunLengthEncodeComplex<U16>;
			this->shiftSize_ = sizeof(U16)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width_ = sizeof(U16);
		} else if(this->VCFheader_.size() <= Constants::UPPER_LIMIT_SAMPLES_32B - 1){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " < " << Constants::UPPER_LIMIT_SAMPLES_32B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 32-bit width..." << std::endl;
			this->encode_ = &TomahawkImportRLE::RunLengthEncodeSimple<U32>;
			this->encodeComplex_ = &TomahawkImportRLE::RunLengthEncodeComplex<U32>;
			this->shiftSize_ = sizeof(U32)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH;
			this->bit_width_ = sizeof(U32);
		} else {
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " > " << Constants::UPPER_LIMIT_SAMPLES_32B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->VCFheader_.size() << " < " << Constants::UPPER_LIMIT_SAMPLES_64B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 64-bit width..." << std::endl;
			this->encode_ = &TomahawkImportRLE::RunLengthEncodeSimple<U64>;
			this->encodeComplex_ = &TomahawkImportRLE::RunLengthEncodeComplex<U64>;
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

	inline const BYTE& getBitWidth(void) const{ return this->bit_width_; }

private:
	template <class T> void RunLengthEncodeSimple (const VCF::VCFLine& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs);
	template <class T> void RunLengthEncodeComplex(const VCF::VCFLine& line, IO::BasicBuffer& meta, IO::BasicBuffer& runs);

private:
	VCF::VCFHeader& VCFheader_;
	rleFunction encode_;			// encoding function
	rleFunction encodeComplex_;		// encoding function
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
	for(U32 i = 1; i < this->VCFheader_.size(); ++i){
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


			run_length = 1;
			type = curType;
			++runsCount;

		} else ++run_length;
	}

	// Encode final
	__dump =  (run_length & (((T)1 << this->shiftSize_) - 1)) << Constants::TOMAHAWK_SNP_PACK_WIDTH;
	__dump ^= (curType & ((1 << Constants::TOMAHAWK_SNP_PACK_WIDTH) - 1));
	runs += __dump;
	total_samples += run_length;
	//runs.pointer += PACK3(curType, &runs.data[runs.pointer], run_length);
	this->helper_[type] += run_length;
	++runsCount;

	if(total_samples != this->VCFheader_.size()){
		std::cerr << Helpers::timestamp("ERROR", "RLE") << "Sum of run lengths does not equal number of samples: " << total_samples << "/" << this->VCFheader_.size() << std::endl;
		exit(1);
	}

	this->helper_.calculateMAF();
	this->helper_.calculateHardyWeinberg();

	// Position
	U32& position = *reinterpret_cast<U32*>(&meta[meta.pointer - 5]);
	position <<= 2;
	position |= this->helper_.phased << 1;
	position |= this->helper_.missingValues << 0;
	meta += this->helper_.MAF;
	meta += this->helper_.HWE_P;

	/*
	if(runs.pointer - startOffset > VCFheader_.size()/2){
		this->savings += (runs.pointer - startOffset) - VCFheader_.size()/2;
		std::cerr << runs.pointer - startOffset << '\t' << (runs.pointer - startOffset) - VCFheader_.size()/2 << '\t' << VCFheader_.size()/2 << '\t' << this->helper_.MAF << std::endl;
		runs.pointer = startOffset;

		U32 i = 0;
		for(; i < line.simple_.size(); i += 2){
			BYTE add = (Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[i].snpA - 46] << 6) |
					   (Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[i].snpB - 46] << 4) |
					   (Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[i+1].snpA - 46] << 2) |
					   (Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[i+1].snpB - 46]);
			runs += add;
		}

		BYTE add = 0;
		for(; i < line.simple_.size(); ++i){
			add <<= 2;
			add |= (Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[i].snpA - 46] << 2) | (Constants::TOMAHAWK_ALLELE_LOOKUP[line.simple_[i].snpB - 46]);
		}
		runs += add;
		runsCount = i;
	}
	*/
	meta += runsCount;

	//std::cerr << this->helper_.MAF << '\t' << runsCount << std::endl;

	this->helper_.reset();
	//std::cerr << "Sum: " << sum << '/' << this->VCFheader_.size() << std::endl;
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
	for(U32 i = 1; i < this->VCFheader_.size(); ++i){
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

	if(total_samples != this->VCFheader_.size()){
		std::cerr << Helpers::timestamp("ERROR", "RLE") << "Sum of run lengths does not equal number of samples: " << total_samples << "/" << this->VCFheader_.size() << std::endl;
		exit(1);
	}

	this->helper_.calculateMAF();
	this->helper_.calculateHardyWeinberg();

	// Position
	U32& position = *reinterpret_cast<U32*>(&meta[meta.pointer - 5]);
	//std::cerr << "Position start " << position << std::endl;
	position <<= 2;
	position |= this->helper_.phased << 1;
	position |= this->helper_.missingValues << 0;
	//std::cerr << "Position end " << position << std::endl;


	meta += this->helper_.MAF;
	meta += this->helper_.HWE_P;
	meta += runsCount;

	this->helper_.reset();
	//std::cerr << "Sum: " << sum << '/' << this->VCFheader_.size() << std::endl;
}

}
}

#endif /* TOMAHAWKIMPORTRLE_H_ */
