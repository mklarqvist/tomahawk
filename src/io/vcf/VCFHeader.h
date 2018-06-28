#ifndef VCFHEADER_H_
#define VCFHEADER_H_

#include <algorithm>

#include "support/helpers.h"
#include "VCFHeaderConstants.h"
#include "VCFHeaderContig.h"
#include "VCFHeaderLine.h"
#include "io/reader.h"
#include "algorithm/open_hashtable.h"
#include "io/basic_buffer.h"
#include "io/compression/tgzf_controller.h"

namespace Tomahawk {
namespace VCF{

class VCFHeader {
	typedef VCFHeader self_type;
	typedef Tomahawk::Hash::HashTable<std::string, S32> hash_table;
	typedef VCFHeaderContig contig_type;
	typedef IO::TGZFController tgzf_type;
	typedef IO::BasicBuffer buffer_type;

	enum VCF_ERROR_TYPE {VCF_PASS, VCF_ERROR_LINE1, VCF_ERROR_LINES, VCF_ERROR_SAMPLE, STREAM_BAD};

public:
	VCFHeader();
	~VCFHeader();

	void unsetBorrowedPointers(void){
		this->contigsHashTable = nullptr;
		this->sampleHashTable = nullptr;
	}

	inline bool good(void) const{ return(this->error_bit == VCF_PASS); }
	inline bool valid(void) const{ return(this->version > 0); }
	inline void setVersion(float version){ this->version = version; }
	inline const float& getVersion(void) const{ return(this->version); }
	inline U32 getContigs(void) const{ return this->contigs.size(); }
	inline const contig_type& operator[](const U32 p) const{ return(this->contigs[p]); }
	inline contig_type& getContig(const U32 p){ return this->contigs[p]; }
	inline U32 getLines(void) const{ return this->lines.size(); }
	inline const U64& size(void) const{ return this->samples; }

	inline bool getContig(const std::string& contig, S32*& retValue) const{
		return(this->contigsHashTable->GetItem(&contig[0], &contig, retValue, contig.size()));
	}

	inline bool getSample(const std::string& sample, S32*& retValue) const{
		return(this->sampleHashTable->GetItem(&sample[0], &sample, retValue, sample.size()));
	}

	bool parse(reader& stream);
	bool parse(const char* const data, const U32& length);

	bool writeTGZFLiterals(std::ofstream& stream) const{
		buffer_type temp;
		tgzf_type tgzf_controller;

		for(U32 i = 0; i < this->literal_lines.size(); ++i){
			std::cerr << this->literal_lines[i] << std::endl;
			temp += this->literal_lines[i];
		}
		tgzf_controller.Deflate(temp);
		stream.write(temp.data(), temp.size());
		tgzf_controller.Clear();
		temp.deleteAll();

		return true;
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
	void buildSampleTable(U64 samples);
	bool __parseFirstLine(reader& stream);
	bool __parseHeaderLines(reader& stream);
	bool __parseSampleLine(reader& stream);

	bool __parseFirstLine(const char* const data, U32& offset);
	bool __parseHeaderLines(const char* const data, U32& offset);
	bool __parseSampleLine(const char* const data, U32& offset, const U32& length);

public:
	VCF_ERROR_TYPE error_bit;				// parse error bit
	U64 samples;							// number of samples
	float version;							// VCF version
	std::string literal;					// string copy of header data
	std::vector<contig_type> contigs;		// contigs
	std::vector<std::string> sampleNames;	// sample names
	std::vector<VCFHeaderLine> lines;		// header lines
	std::vector<std::string> literal_lines; // vcf line literals
	hash_table* contigsHashTable;			// hash table for contig names
	hash_table* sampleHashTable;			// hash table for sample names
};

}
} /* namespace Tomahawk */

#endif /* VCFHEADER_H_ */
