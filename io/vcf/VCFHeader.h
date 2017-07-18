#ifndef VCFHEADER_H_
#define VCFHEADER_H_

#include <algorithm>

#include "../../helpers.h"
#include "VCFHeaderConstants.h"
#include "../../algorithm/OpenHashTable.h"

namespace Tomahawk {
namespace VCF{

struct VCFHeaderContig{
	typedef VCFHeaderContig self_type;

public:
	VCFHeaderContig() : length(0), tomahawkBlocks(0){}
	~VCFHeaderContig(){}

	friend std::ostream& operator<<(std::ostream& out, const self_type& contig){
		out << contig.name << '\t' << contig.length;
		return(out);
	}

	void operator++(void){ ++this->tomahawkBlocks; }

public:
	std::string name;
	U32 length;
	// keep track of how many blocks we've seen for this contig
	// used during import
	U32 tomahawkBlocks;
};

struct VCFHeaderLineKeyValue{
	typedef VCFHeaderLineKeyValue self_type;

public:
	VCFHeaderLineKeyValue(){}
	~VCFHeaderLineKeyValue(){}

	friend std::ostream& operator<<(std::ostream& out, const self_type& pair){
		out << pair.lKEY << '/' << pair.lVALUE << '\t' << std::string(&pair.KEY[0], pair.lKEY) << '\t' << std::string(&pair.VALUE[0], pair.lVALUE);
		return(out);
	}

public:
	U32 lKEY;
	U32 lVALUE;
	const char* KEY;
	const char* VALUE;
};

struct VCFHeaderLine{
public:
	VCFHeaderLine(const char* data, const U32 size) : size_(size), data_(data){}
	~VCFHeaderLine(){}

	const U32 size(void) const{ return this->pairs_.size(); }
	const VCFHeaderLineKeyValue& operator[](const U32 p) const{ return this->pairs_[p]; }

	bool isValid(void) const{ return(this->size_ > 2 && (this->data_[0] == '#' && this->data_[1] == '#')); }
	bool isCONTIG(void) const{ return(strncasecmp(&Constants::HEADER_CONTIG[0], this->data_, Constants::HEADER_CONTIG.size()) == 0); }

	friend std::ostream& operator<<(std::ostream& out, const VCFHeaderLine& pair){
		for(U32 i = 0; i < pair.pairs_.size(); ++i)
			out << pair[i] << '\n';

		return(out);
	}

	bool Parse(void){
		if(!this->isValid()){
			std::cerr << "Invalid line" << std::endl;
			return false;
		}

		const char* match = std::find(this->data_, &this->data_[this->size_], '=');
		if(*match != '='){
			std::cerr << "invalid entry no equal match" << std::endl;
			return false;
		}
		U32 matchPos = match - this->data_ + 1;
		if(this->data_[matchPos] == '<'){
			if(this->data_[this->size_] != '>'){
				std::cerr << "Illegal header line: " << this->data_[this->size_] << std::endl;
				return false;
			}

			++matchPos;
			while(this->nextKey(matchPos)){

			}
		} else {
			//Todo: this value is just text
			//std::cerr << "value is just text" << std::endl;
		}

		return true;
	}

private:
	bool nextKey(U32& startPos){
		if(this->data_[startPos] == '>')
			return false;

		//std::cerr << "Searching from " << startPos << " / " << this->size_ << std::endl;

		const char* match = std::find(&this->data_[startPos], &this->data_[this->size_], '=');
		if(*match != '='){
			std::cerr << "invalid entry no equal match in next key" << std::endl;
			return false;
		}
		U32 matchPos = match - this->data_;
		VCFHeaderLineKeyValue entry;
		entry.KEY = &this->data_[startPos];
		entry.lKEY = matchPos - startPos;
		//std::cerr << "Startpos: " << startPos << "->" << matchPos - 1 << " length " << matchPos - startPos - 1 << std::endl;

		startPos = matchPos + 1;

		char match_token = ',';
		BYTE adjust_value = 0;
		if(this->data_[startPos] == '"'){
			//std::cerr << "Search for quotes at pos " << startPos << std::endl;
			match_token = '"';
			adjust_value = 1;
		}

		match = std::find(&this->data_[startPos + adjust_value], &this->data_[this->size_], match_token);
		//std::cerr << "match is " << *match << " at " << match - this->data_ << std::endl;
		if(*match == '>'){
			entry.VALUE = &this->data_[startPos];
			entry.lVALUE = this->size_ - startPos;
			startPos = matchPos + 1;
			this->pairs_.push_back(entry);
			return false;
		} else if(*match != match_token){
			std::cerr << "invalid entry no comma match" << std::endl;
			return false;
		}
		matchPos = match - this->data_;
		entry.VALUE = &this->data_[startPos];
		entry.lVALUE = matchPos - startPos + adjust_value;
		startPos = matchPos + 1;
		this->pairs_.push_back(entry);
		return true;
	}

private:
	const U32 size_;
	const char* data_;
	std::vector<VCFHeaderLineKeyValue> pairs_;
};

class VCFHeader {
	typedef Tomahawk::Hash::HashTable<std::string, U32> hashtable;

public:
	VCFHeader() : contigsHashTable_(nullptr), sampleHashTable_(nullptr){}
	~VCFHeader(){ delete this->contigsHashTable_; }

	bool valid(void) const{ return(this->version_ > 0); }
	void setVersion(S32 version){ this->version_ = version; }
	inline void addContig(const std::string& contig, U32 value){
		this->contigsHashTable_->SetItem(&contig[0], &contig, value, contig.size());
	}

	inline bool getContig(const std::string& contig, U32*& retValue){
		return(this->contigsHashTable_->GetItem(&contig[0], &contig, retValue, contig.size()));
	}

	inline void addSample(const std::string& sample){
		this->sampleNames_.push_back(sample);
		this->sampleHashTable_->SetItem(&sample[0], &sample, this->sampleNames_.size()-1, sample.size());
	}

	inline bool getSample(const std::string& sample, U32*& retValue){
		return(this->sampleHashTable_->GetItem(&sample[0], &sample, retValue, sample.size()));
	}

	U32 getContigs(void) const{ return this->contigs_.size(); }
	VCFHeaderContig& getContig(const U32 p){ return this->contigs_[p]; }
	U32 getLines(void) const{ return this->lines_.size(); }

	bool checkLine(const char* data, const U32 length){
		VCFHeaderLine line(data, length);
		if(line.Parse()){
			//std::cerr << line;
			this->lines_.push_back(line);
			if(line.isCONTIG()){
				VCFHeaderContig contig;
				BYTE found = 0;
				for(U32 i = 0; i < line.size(); ++i){
					if(strncmp(line[i].KEY, "ID", 2) == 0){
						contig.name = std::string(line[i].VALUE, line[i].lVALUE);
						//std::cerr << std::string(line[i].VALUE, line[i].lVALUE) << std::endl;
						++found;
					} else if(strncmp(line[i].KEY, "length", 6) == 0){
						contig.length = atoi(line[i].VALUE);
						//std::cerr << contig.length << std::endl;
						++found;
					}
				}

				if(found != 2){
					std::cerr << Helpers::timestamp("WARNING","VCF") << "Illegal contig entry line with no length defined!" << std::endl;
					std::cerr << Helpers::timestamp("WARNING","VCF") << "Offending line: " << std::string(data, length+1) << std::endl;
					contig.length = ~0;
					//return false;
				}
				this->contigs_.push_back(contig);
				return true;
			}

			return true;
		}

		std::cerr << Helpers::timestamp("ERROR","VCF") << "Failed to parse VCF..." << std::endl;
		return false;
	}

	bool BuildContigTable(void){
		U32* retValue;

		if(this->contigs_.size() < 1024)
			this->contigsHashTable_ = new hashtable(1024);
		else
			this->contigsHashTable_ = new hashtable(this->contigs_.size() * 2);

		std::cerr << Helpers::timestamp("LOG", "VCF") << "Constructing lookup table for " << this->contigs_.size() << " contigs..." << std::endl;


		for(U32 i = 0; i < this->contigs_.size(); ++i){
			//std::cerr << i << '/' << this->contigs_.size() << std::endl;
			if(!(*this).getContig(this->contigs_[i].name, retValue)){
				//std::cerr << "not in here. inserting" << std::endl;
				(*this).addContig(this->contigs_[i].name, i);
			} else {
				std::cerr << Helpers::timestamp("ERROR", "VCF") << "Duplicated contig found (" << this->getContig(*retValue).name << "). Illegal..." << std::endl;
				return false;
			}
		}
		return true;
	}

	inline void setSamples(U32 samples){
		this->samples_ = samples;
		delete this->sampleHashTable_;

		if(this->samples_ < 1024)
			this->sampleHashTable_ = new hashtable(1024);
		else
			this->sampleHashTable_ = new hashtable(this->samples_ * 2);
	}
	inline const U64 size(void) const{ return this->samples_; }

	inline const VCFHeaderContig& operator[](const U32 p) const{ return(this->contigs_[p]); }

public:
	U64 samples_;
	S32 version_;
	std::vector<VCFHeaderLine> lines_;
	std::vector<VCFHeaderContig> contigs_;
	std::vector<std::string> sampleNames_;
	hashtable* contigsHashTable_;
	hashtable* sampleHashTable_;
};

}
} /* namespace Tomahawk */

#endif /* VCFHEADER_H_ */
