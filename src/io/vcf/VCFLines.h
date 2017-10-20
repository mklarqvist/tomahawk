#ifndef VCFLINES_H_
#define VCFLINES_H_

#include <iostream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <vector>
#include <algorithm>

#include "../../support/helpers.h"
#include "../../support/MagicConstants.h"

namespace Tomahawk{
namespace VCF{

//
#pragma pack(1)
struct VCFDiploidGenotype{
public:
	VCFDiploidGenotype(); // Has no ctor or dtor
	~VCFDiploidGenotype();

	bool hasMissing(void) const{return(this->snpA == '.' || this->snpB == ','); }

	char snpA;
	char separator;
	char snpB;
	char spacer;

	friend std::ostream& operator<<(std::ostream& stream, const VCFDiploidGenotype& entry){
		stream << entry.snpA << entry.separator << entry.snpB;
		return stream;
	}
};

class VCFLineDataInterface{
public:
	VCFLineDataInterface(const U32 samples) : samples_(samples), dataLength_(0), inputData_(nullptr){}
	virtual ~VCFLineDataInterface(){}
	virtual bool Parse(void) =0;
	void SetData(const char* data, const U32 dataLength){ this->inputData_ = data; this->dataLength_ = dataLength; }
	//virtual const VCFDiploidGenotype& operator[](const U32 position) const =0;
	const U32& size(void) const{ return this->samples_; }

protected:
	const U32 samples_;
	U32 dataLength_;
	const char* inputData_;
};

class VCFLineDataSimple : public VCFLineDataInterface{
public:
	VCFLineDataSimple(const U32 samples) : VCFLineDataInterface(samples), data_(nullptr){}
	~VCFLineDataSimple(){}

	const VCFDiploidGenotype& operator[](const U32 position) const{ return this->data_[position]; }

	bool Parse(void){ // From interface
		this->data_ = reinterpret_cast<const VCFDiploidGenotype*>(this->inputData_);
		return true;
	}

public:
	const VCFDiploidGenotype* data_;
};

class VCFLineDataComplex : public VCFLineDataInterface {
public:
	VCFLineDataComplex(const U32 samples) : VCFLineDataInterface(samples), data_(samples){}
	~VCFLineDataComplex(){}

	const VCFDiploidGenotype* operator[](const U32 position) const{ return this->data_[position]; }
	bool Parse(void){
		this->data_.clear();
		uint32_t search_position = 0;
		uint32_t delimiters_found = 0;
		while(true){ // while there is samples in line
			const char* found = std::find(&this->inputData_[search_position], &this->inputData_[this->dataLength_], Constants::VCF_DELIMITER);

			if(*found != Tomahawk::VCF::Constants::VCF_DELIMITER){
				//std::cerr << "break no match " << (*found == '\n') << '\t' << (int)*found << '\t' << *found << '\t' << found - this->inputData_ << '/' << this->dataLength_ << std::endl;
				break;
			}

			//std::cerr << std::string(&this->inputData_[search_position], (found - this->inputData_ + 1) - search_position) << std::endl;
			this->data_.push_back(reinterpret_cast<const VCFDiploidGenotype*>(&this->inputData_[search_position]));
			//std::cerr << this->data_[this->data_.size()-1]->snpA << '/' << this->data_[this->data_.size()-1]->snpB << std::endl;
			search_position = found - this->inputData_ + 1;
			++delimiters_found;
		}

		//std::cerr << std::string(&this->inputData_[search_position], this->dataLength_ - search_position + 1) << std::endl;
		//std::cerr << "final: " << this->dataLength_ - search_position + 1 << '\t' << this->dataLength_ << '\t' << search_position << std::endl;
		this->data_.push_back(reinterpret_cast<const VCFDiploidGenotype*>(&this->inputData_[search_position]));
		++delimiters_found;

		if(delimiters_found == this->samples_)
			return true;

		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Found " << delimiters_found << " samples in line but expected " << this->samples_ << "..." << std::endl;
		exit(1);
		return false;

	} // From interface

public:
	std::vector<const VCFDiploidGenotype*> data_;
};

class VCFLine{
public:
	VCFLine(const U32 samples): simple_(samples), complex_(samples){}
	~VCFLine(void){}

	inline const bool checkSeparator(const char& separator) const{
		if(separator == '/')
			return true;
		else if(separator == '|')
			return true;
		else return false;
	}


	bool Parse(const char* source, const U32 sourceLength){
		BYTE found = 0;
		U32 sourceLastPosition = 0;
		U32 sourceFoundPosition = 0;

		while(true){
			const char* match = std::find(&source[sourceLastPosition], &source[sourceLength], Tomahawk::VCF::Constants::VCF_DELIMITER);
			if(*match != Tomahawk::VCF::Constants::VCF_DELIMITER){
				std::cerr << Helpers::timestamp("ERROR", "VCF") << "Illegal VCF line" << std::endl;
				return false;
			}

			sourceFoundPosition = match - source;

			switch(found){
			case 0:
				this->CHROM = &source[0];
				this->lCHROM = sourceFoundPosition - sourceLastPosition;
				break;
			case 1:
				this->POS = &source[sourceLastPosition];
				this->lPOS = sourceFoundPosition - sourceLastPosition;
				this->position = atoi(this->POS);
				break;
			case 2:
				this->ID = &source[sourceLastPosition];
				this->lID = sourceFoundPosition - sourceLastPosition;
				break;
			case 3:
				this->REF = &source[sourceLastPosition];
				this->lREF = sourceFoundPosition - sourceLastPosition;
				break;
			case 4:
				this->ALT = &source[sourceLastPosition];
				this->lALT = sourceFoundPosition - sourceLastPosition;
				break;
			case 5:
				this->QUAL = &source[sourceLastPosition];
				this->lQUAL = sourceFoundPosition - sourceLastPosition;
				break;
			case 6:
				this->FILTER = &source[sourceLastPosition];
				this->lFILTER = sourceFoundPosition - sourceLastPosition;
				break;
			case 7:
				this->INFO = &source[sourceLastPosition];
				this->lINFO = sourceFoundPosition - sourceLastPosition;
				break;
			case 8:
				this->FORMAT = &source[sourceLastPosition];
				this->lFORMAT = sourceFoundPosition - sourceLastPosition;
				break;
			}

			sourceLastPosition = sourceFoundPosition + 1;
			++found;

			if(found == 9){
				this->isComplex();
				this->SetReference();
				if(!this->getComplex() && this->IsSimple()){
					this->simple_.SetData(&source[sourceLastPosition], sourceLength - sourceLastPosition - 1);
					this->simple_.Parse();
				} else if(this->getComplex() && this->IsSimple()){
					this->complex_.SetData(&source[sourceLastPosition], sourceLength - sourceLastPosition - 1);
					this->complex_.Parse();
				}
				return true;
			}
		}

		return false;
	}

	// Complex: defined as FORMAT field equals "GT" && site is biallelic
	const bool& getComplex(void) const{ return this->Complex; }

	template <class T>
	const float getMissingness(const T& samples) const{
		U64 total = 0;
		if(this->getComplex()){
			for(U32 i = 0; i < samples; ++i){
				if(!this->checkSeparator(this->complex_.data_[i]->separator))
					return(2);

				if(this->complex_.data_[i]->snpA == '.' || this->complex_.data_[i]->snpB == '.')
					++total;
			}
		} else {
			for(U32 i = 0; i < samples; ++i){
				if(!this->checkSeparator(this->simple_.data_[i].separator))
					return(2);

				if(this->simple_[i].snpA == '.' || this->simple_.data_[i].snpB == '.')
					++total;
			}
		}
		return((float)total/samples);
	}

	bool isComplex(void){
		if(strncmp(this->FORMAT, &Tomahawk::VCF::Constants::GT_ONLY[0], Tomahawk::VCF::Constants::GT_ONLY.size()) == 0 && this->lFORMAT == 2)
			this->Complex = false;
		else {
			if(strncmp(this->FORMAT, &Tomahawk::VCF::Constants::GT_ONLY[0], Tomahawk::VCF::Constants::GT_ONLY.size()) == 0)
				this->Complex = true;
			else {
				std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not parse GT information..." << std::endl;
				exit(1);
			}
		}

		return this->Complex;
	}
	const bool IsSimple(void) const{ return(this->lALT == 1 && this->lREF == 1); }

	void SetReference(void){
		this->ref_alt = 0;

		switch(this->REF[0]){
		case 'A': this->ref_alt ^= Tomahawk::Constants::REF_ALT_A << 4; break;
		case 'T': this->ref_alt ^= Tomahawk::Constants::REF_ALT_T << 4; break;
		case 'G': this->ref_alt ^= Tomahawk::Constants::REF_ALT_G << 4; break;
		case 'C': this->ref_alt ^= Tomahawk::Constants::REF_ALT_C << 4; break;
		case '.': this->ref_alt ^= Tomahawk::Constants::REF_ALT_N << 4; break;
		}

		switch(this->ALT[0]){
		case 'A': this->ref_alt ^= Tomahawk::Constants::REF_ALT_A << 0; break;
		case 'T': this->ref_alt ^= Tomahawk::Constants::REF_ALT_T << 0; break;
		case 'G': this->ref_alt ^= Tomahawk::Constants::REF_ALT_G << 0; break;
		case 'C': this->ref_alt ^= Tomahawk::Constants::REF_ALT_C << 0; break;
		case '.': this->ref_alt ^= Tomahawk::Constants::REF_ALT_N << 0; break;
		}
	}

	const BYTE& getReference(void) const{ return this->ref_alt; }

public:
	U16 lCHROM;
	U16 lPOS;
	U16 lID;
	U16 lREF;
	U16 lALT;
	U16 lQUAL;
	U16 lFILTER;
	U16 lINFO;
	U16 lFORMAT;
	U32 position;
	bool Complex;
	const char* CHROM;
	const char* POS;
	const char* ID;
	const char* REF;
	const char* ALT;
	const char* QUAL;
	const char* FILTER;
	const char* INFO;
	const char* FORMAT;
	BYTE ref_alt;
	VCFLineDataSimple simple_;
	VCFLineDataComplex complex_;
};

}
}



#endif /* VCFLINES_H_ */
