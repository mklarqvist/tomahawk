#ifndef VCF_VCFHEADERLINE_H_
#define VCF_VCFHEADERLINE_H_

#include <algorithm> // for std::find

namespace Tomahawk{
namespace VCF{

struct VCFHeaderLine{
private:
	typedef VCFHeaderLine self_type;

	// Internal helper struct
	struct VCFHeaderLineKeyValue{
		typedef VCFHeaderLineKeyValue self_type;

		VCFHeaderLineKeyValue(){}
		~VCFHeaderLineKeyValue(){}

		friend std::ostream& operator<<(std::ostream& out, const self_type& pair){
			out << pair.KEY << '\t' << pair.VALUE;
			return(out);
		}

		std::string KEY;
		std::string VALUE;
	};
	typedef VCFHeaderLineKeyValue key_value;

public:
	VCFHeaderLine(const char* data, const U32 size) : size_(size), data_(data){}
	~VCFHeaderLine(){}

	inline const U32 size(void) const{ return this->pairs_.size(); }
	inline const key_value& operator[](const U32 p) const{ return this->pairs_[p]; }
	inline bool isValid(void) const{ return(this->size_ > 2 && (this->data_[0] == '#' && this->data_[1] == '#')); }
	inline bool isCONTIG(void) const{
		return(strncasecmp(&Constants::HEADER_CONTIG[0], &this->data_[0], Constants::HEADER_CONTIG.size()) == 0);
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& pair){
		out << pair.data_ << '\n';
		for(U32 i = 0; i < pair.pairs_.size(); ++i)
			out << i << '/' << pair.size() << '\t' << pair[i] << '\n';

		return(out);
	}

	bool Parse(void){
		// Make sure this is a valid VCF header line
		// Rule: has to start with ##
		if(!this->isValid()){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Invalid VCF header line..." << std::endl;
			return false;
		}

		// Attempt to find an equal sign
		const char* match = std::find(this->data_, &this->data_[this->size_], '=');
		if(*match != '='){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: no equal match..." << std::endl;
			return false;
		}

		U32 matchPos = match - this->data_ + 1;
		if(this->data_[matchPos] == '<'){
			if(this->data_[this->size_] != '>'){
				std::cerr << Helpers::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: " << this->data_[this->size_] << std::endl;
				return false;
			}

			++matchPos;
			// Sweep over and assert it is valid
			while(this->nextKey(matchPos)){
				// nothing in body
			}

			// Todo
			//for(U32 i = 0; i < this->pairs_.size(); ++i){
			//	std::cerr << i << '\t' << this->pairs_[i].KEY << '\t' << this->pairs_[i].VALUE << std::endl;
			//}

		} else {
			//Todo: this value is just text
		}
		return true;
	}

private:
	bool nextKey(U32& startPos){
		if(this->data_[startPos] == '>')
			return false;

		const char* match = std::find(&this->data_[startPos], &this->data_[this->size_], '=');
		if(*match != '='){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: no equal match in next key..." << std::endl;
			return false;
		}
		U32 matchPos = match - this->data_;
		VCFHeaderLineKeyValue entry;
		entry.KEY = std::string(&this->data_[startPos], matchPos - startPos);
		//entry.lKEY = matchPos - startPos;
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
			entry.VALUE = std::string(&this->data_[startPos],this->size_ - startPos);
			//entry.lVALUE = this->size_ - startPos;
			startPos = matchPos + 1;
			this->pairs_.push_back(entry);
			return false;
		} else if(*match != match_token){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: no comma match in next key..." << std::endl;
			return false;
		}

		matchPos = match - this->data_;
		entry.VALUE = std::string(&this->data_[startPos], matchPos - startPos + adjust_value);
		//entry.lVALUE = matchPos - startPos + adjust_value;
		startPos = matchPos + 1;
		this->pairs_.push_back(entry);
		return true;
	}

public:
	U32 size_;
	const char* data_;
	std::vector<key_value> pairs_;
};

}
}

#endif /* VCF_VCFHEADERLINE_H_ */
