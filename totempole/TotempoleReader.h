#ifndef TOMAHAWK_TOTEMPOLEREADER_H_
#define TOMAHAWK_TOTEMPOLEREADER_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "../tomahawk/MagicConstants.h"
#include "../TypeDefinitions.h"
#include "../helpers.h"
#include "TotempoleEntry.h"
#include "../algorithm/OpenHashTable.h"
#include "../io/totempole/TotempoleContig.h"
#include "../io/totempole/TotempoleHeader.h"

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
	typedef TotempoleEntry  entry_type;

	typedef Tomahawk::Hash::HashTable<std::string, U32> hashtable;

public:
	TotempoleReader();
	~TotempoleReader();
	bool Open(const std::string filename);
	std::vector<U32> findOverlaps(const Interval& interval) const;

	inline const U32& getLargestBlockSize(void) const{ return this->header.largest_uncompressed; }
	inline const U32& getBlocks(void) const{ return this->header.blocks; }
	inline const U64& getSamples(void) const{ return this->header.samples; }
	inline const U32& getContigs(void) const{ return this->n_contigs; }
	inline const U32& size(void) const{ return this->n_contigs; }
	inline const entry_type& operator[](const U32 p) const{ return this->entries[p]; }
	inline const contig_type& getContig(const U32 contigID) const{ return this->contigs[contigID]; }
	inline const header_type& getHeader(void) const{ return(this->header); }

private:
	bool Validate(std::ifstream& in) const;
	bool BuildHashTables(void);

private:
	std::string filename;	// filename
	U32 filesize;			// filesize
	U32 n_contigs;			// number of contigs
	header_type header;		// header information
	contig_type* contigs;	// contig data
	std::string* samples;	// sample names
	entry_type* entries;	// totempole entries data
	hashtable* contigsHashTable;	// contig name hash table
	hashtable* sampleHashTable;		// smaple name hash table
};

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOTEMPOLEREADER_H_ */
