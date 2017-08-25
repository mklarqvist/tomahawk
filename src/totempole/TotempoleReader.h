#ifndef TOMAHAWK_TOTEMPOLEREADER_H_
#define TOMAHAWK_TOTEMPOLEREADER_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "../support/MagicConstants.h"
#include "../support/TypeDefinitions.h"
#include "../support/helpers.h"
#include "TotempoleEntry.h"
#include "../algorithm/OpenHashTable.h"
#include "../totempole/TotempoleContig.h"
#include "../totempole/TotempoleHeader.h"
#include "../io/TGZFController.h"

namespace Tomahawk {

// Todo: temp interval
struct Interval{
	Interval(U32 contig, U32 from, U32 to) : contigID(contig), from(from), to(to){}
	~Interval(){}

	U32 contigID;
	U32 from;
	U32 to;
};

class TotempoleReader {
	typedef Totempole::TotempoleHeader header_type;
	typedef Totempole::TotempoleContigBase contig_base_type;
	typedef Totempole::TotempoleContig contig_type;
	typedef TotempoleEntry entry_type;
	typedef Tomahawk::Hash::HashTable<std::string, S32> hash_table;
	typedef IO::TGZFHeader tgzf_type;
	typedef IO::TGZFController tgzf_controller_type;
	typedef IO::BasicBuffer buffer_type;

public:
	TotempoleReader();
	~TotempoleReader();
	bool Open(const std::string filename);
	std::vector<U32> findOverlaps(const Interval& interval) const;
	const entry_type& front(void) const{ return(this->entries[0]); }
	const entry_type& back(void) const{ return(this->entries[this->getBlocks() - 1]); };

	inline const U32& getLargestBlockSize(void) const{ return this->header.largest_uncompressed; }
	inline const U32& getBlocks(void) const{ return this->header.blocks; }
	inline const U64& getSamples(void) const{ return this->header.samples; }
	inline const U32& getContigs(void) const{ return this->n_contigs; }
	inline const U32& size(void) const{ return this->n_contigs; }
	inline const entry_type& operator[](const U32 p) const{ return this->entries[p]; }
	inline const contig_type& getContig(const U32 contigID) const{ return this->contigs[contigID]; }
	inline bool getContig(const std::string& string, S32*& contigID) const{
		if(this->contigsHashTable->GetItem(&string[0], &string, contigID, string.size()))
			return true;

		return false;
	}
	inline const header_type& getHeader(void) const{ return(this->header); }
	inline const contig_base_type* getContigBase(const U32 contigID) const{ return(reinterpret_cast<const contig_base_type*>(&this->contigs[contigID])); }

	hash_table* getContigHTablePointer(void) const{ return(this->contigsHashTable); }
	hash_table* getSampleHTablePointer(void) const{ return(this->sampleHashTable); }
	bool writeLiterals(std::ofstream& stream);

private:
	bool Validate(std::ifstream& in) const;
	void BuildUpdateContigs(void);
	bool ValidateEOF(std::ifstream& in);
	bool BuildHashTables(void);

public:
	std::ifstream stream;	// filestream
	std::string filename;	// filename
	U32 filesize;			// filesize
	U32 n_contigs;			// number of contigs
	header_type header;		// header information
	std::string literals;	// literal data
	tgzf_controller_type tgzf_controller;	// tgzf controller
	contig_type* contigs;	// contig data
	std::string* samples;	// sample names
	entry_type* entries;	// totempole entries data
	hash_table* contigsHashTable;	// contig name hash table
	hash_table* sampleHashTable;	// sample name hash table
};

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOTEMPOLEREADER_H_ */
