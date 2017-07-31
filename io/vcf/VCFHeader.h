#ifndef VCFHEADER_H_
#define VCFHEADER_H_

#include <algorithm>

#include "../../helpers.h"
#include "VCFHeaderConstants.h"
#include "../reader.h"
#include "../../algorithm/OpenHashTable.h"

namespace Tomahawk {
namespace VCF{

struct VCFHeaderContig{
	typedef VCFHeaderContig self_type;

public:
	VCFHeaderContig() : length(0), tomahawkBlocks(0){}
	~VCFHeaderContig(){}

	inline void operator++(void){ ++this->tomahawkBlocks; }
	inline void operator--(void){ --this->tomahawkBlocks; }
	template <class T> inline void operator+=(const T value){ this->tomahawkBlocks += value; }
	template <class T> inline void operator-=(const T value){ this->tomahawkBlocks -= value; }

	friend std::ostream& operator<<(std::ostream& out, const self_type& contig){
		out << contig.name << '\t' << contig.length;
		return(out);
	}

public:
	std::string name;
	U32 length;
	// keep track of how many blocks we've seen for this contig
	// used during import
	U32 tomahawkBlocks;
};

struct VCFHeaderLine{
private:
	typedef VCFHeaderLine self_type;

	// Internal helper struct
	struct VCFHeaderLineKeyValue{
		typedef VCFHeaderLineKeyValue self_type;

		VCFHeaderLineKeyValue(){}
		~VCFHeaderLineKeyValue(){}

		friend std::ostream& operator<<(std::ostream& out, const self_type& pair){
			out << pair.lKEY << '/' << pair.lVALUE << '\t' << std::string(&pair.KEY[0], pair.lKEY) << '\t' << std::string(&pair.VALUE[0], pair.lVALUE);
			return(out);
		}

		U32 lKEY;
		U32 lVALUE;
		const char* KEY;
		const char* VALUE;
	};
	typedef VCFHeaderLineKeyValue key_value;

public:
	VCFHeaderLine(const char* data, const U32 size) : size_(size), data_(data){}
	~VCFHeaderLine(){}

	inline const U32 size(void) const{ return this->pairs_.size(); }
	inline const key_value& operator[](const U32 p) const{ return this->pairs_[p]; }
	inline bool isValid(void) const{ return(this->size_ > 2 && (this->data_[0] == '#' && this->data_[1] == '#')); }
	inline bool isCONTIG(void) const{ return(strncasecmp(&Constants::HEADER_CONTIG[0], this->data_, Constants::HEADER_CONTIG.size()) == 0); }

	friend std::ostream& operator<<(std::ostream& out, const self_type& pair){
		for(U32 i = 0; i < pair.pairs_.size(); ++i)
			out << pair[i] << '\n';

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
				std::cerr << Helpers::timestamp("ERROR", "VCF") << "Corrupted VCF header entry..." << this->data_[this->size_] << std::endl;
				return false;
			}

			++matchPos;
			// Sweep over and assert it is valid
			while(this->nextKey(matchPos)){
				// nothing in body
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
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: no equal match in next key..." << std::endl;
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
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: no comma match in next key..." << std::endl;
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
	std::vector<key_value> pairs_;
};

class VCFHeader {
	typedef VCFHeader self_type;
	typedef Tomahawk::Hash::HashTable<std::string, U32> hash_table;
	typedef VCFHeaderContig contig_type;

	enum VCF_ERROR_TYPE {VCF_PASS, VCF_ERROR_LINE1, VCF_ERROR_LINES, VCF_ERROR_SAMPLE, STREAM_BAD};

public:
	VCFHeader() : error_bit(VCF_PASS), samples(0), version(0), contigsHashTable(nullptr), sampleHashTable(nullptr){}
	~VCFHeader(){ delete this->contigsHashTable; }

	inline bool good(void) const{ return(this->error_bit == VCF_PASS); }
	inline bool valid(void) const{ return(this->version > 0); }
	inline void setVersion(float version){ this->version = version; }
	inline const float& getVersion(void) const{ return(this->version); }
	inline U32 getContigs(void) const{ return this->contigs.size(); }
	inline contig_type& getContig(const U32 p){ return this->contigs[p]; }
	inline U32 getLines(void) const{ return this->lines.size(); }
	inline const U64& size(void) const{ return this->samples; }
	inline const contig_type& operator[](const U32 p) const{ return(this->contigs[p]); }

	inline bool getContig(const std::string& contig, U32*& retValue) const{
		return(this->contigsHashTable->GetItem(&contig[0], &contig, retValue, contig.size()));
	}

	inline bool getSample(const std::string& sample, U32*& retValue) const{
		return(this->sampleHashTable->GetItem(&sample[0], &sample, retValue, sample.size()));
	}

	friend reader& operator>>(reader& stream, self_type& base){
		if(!base.__parseFirstLine(stream))
			return stream;

		// Read remainder lines
		if(!base.__parseHeaderLines(stream))
			return stream;

		if(!base.buildContigTable())
			return stream;

		// Read samples line
		if(!base.__parseSampleLine(stream))
			return stream;

		return(stream);
	}

private:
	// These functions are unsafe as they require contigHashTable to be
	// set prior to calling
	// no tests are made to check
	inline void addContig(const std::string& contig, U32 value){
		this->contigsHashTable->SetItem(&contig[0], &contig, value, contig.size());
	}

	inline void addSample(const std::string& sample){
		this->sampleNames.push_back(sample);
		this->sampleHashTable->SetItem(&sample[0], &sample, this->sampleNames.size()-1, sample.size());
	}

	// Internal overload helpers
	bool checkLine(const char* data, const U32 length);
	bool buildContigTable(void);
	void buildSampleTable(U32 samples);
	bool __parseFirstLine(reader& stream);
	bool __parseHeaderLines(reader& stream);
	bool __parseSampleLine(reader& stream);

public:
	VCF_ERROR_TYPE error_bit;				// parse error bit
	U64 samples;							// number of samples
	float version;							// VCF version
	std::vector<VCFHeaderLine> lines;		// header lines
	std::vector<contig_type> contigs;		// contigs
	std::vector<std::string> sampleNames;	// sample names
	hash_table* contigsHashTable;			// hash table for contig names
	hash_table* sampleHashTable;			// hash table for sample names
};

}
} /* namespace Tomahawk */

#endif /* VCFHEADER_H_ */
