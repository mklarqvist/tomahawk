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
	TotempoleReader() : filesize(0), n_contigs(0), contigs(nullptr), samples(nullptr), entries(nullptr), contigsHashTable(nullptr), sampleHashTable(nullptr){}
	~TotempoleReader(){
		delete [] this->contigs;
		delete [] this->entries;
		delete [] this->samples;

		delete this->contigsHashTable;
		delete this->sampleHashTable;
	}

	bool Open(const std::string filename){
		if(filename.size() == 0){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "IO") << "No input filename..." << std::endl;
			return false;
		}

		this->filename = filename;

		std::ifstream reader(this->filename, std::ios::in | std::ios::binary | std::ios::ate);
		if(!reader.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "IO") << "Could not open: " << this->filename << "..." << std::endl;
			return false;
		}

		this->filesize = reader.tellg();
		reader.seekg(0);

		if(this->filesize <= 0){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "IO") << "File size is 0..." << std::endl;
			return false;
		}


		if(this->filesize < Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TOTEMPOLE") << "Failed MAGIC..." << std::endl;
			return false;
		}

		// Reader header and validate
		if(!this->Validate(reader)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TOTEMPOLE") << "Could not validate Totempole header..." << std::endl;
			return false;
		}

		// Load header data
		reader >> this->header;
#if DEBUG == 1
		std::cerr << this->header << std::endl;
#endif


		// Get number of contigs
		reader.read(reinterpret_cast<char *>(&this->n_contigs), sizeof(U32));
#if DEBUG == 1
		std::cerr << this->n_contigs << std::endl;
#endif
		this->contigs = new contig_type[this->size()];
		contig_base_type* contig_base = reinterpret_cast<contig_base_type*>(this->contigs);
		for(U32 i = 0; i < this->size(); ++i){
			reader >> contig_base[i];
#if DEBUG == 1
			std::cerr << contig_base[i] << std::endl;
#endif
		}

		exit(1);

		char temp_buffer[65536];
		this->samples = new std::string[this->getSamples()];
		for(U32 i = 0; i < this->getSamples(); ++i){
			reader.read(&temp_buffer[0], sizeof(U32));
			const U32 length = *reinterpret_cast<const U32*>(&temp_buffer[0]);
			reader.read(&temp_buffer[sizeof(U32)], length);
			this->samples[i] = std::string(&temp_buffer[sizeof(U32)], length);
#if DEBUG == 1
			std::cerr << i << '\t' << samples[i] << std::endl;
#endif
		}

		if(reader.tellg() != this->header.offset){
			std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Corrupt file" << std::endl;
			std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << reader.tellg() << '/' << this->header.offset << std::endl;
			return false;
		}

		// Populate Totempole entries
		this->entries = new TotempoleEntry[this->getBlocks()];
		for(U32 i = 0; i < this->getBlocks(); ++i){
			reader >> this->entries[i];
#if DEBUG == 1
			std::cerr << this->entries[i] << std::endl;
#endif
		}

		// Find boundaries for Totempole blocks
		// Master index of indices
		U32 lastContigID = this->entries[0].contigID;
		contigs[lastContigID].minPosition = this->entries[0].minPosition;
		contigs[lastContigID].blocksStart = 0;
		for(U32 i = 1; i < this->getSamples(); ++i){
			if(lastContigID != this->entries[i].contigID){
				contigs[lastContigID].maxPosition = this->entries[i-1].maxPosition;
				contigs[lastContigID].blocksEnd = i;
				contigs[this->entries[i].contigID].minPosition = this->entries[i].minPosition;
				contigs[this->entries[i].contigID].blocksStart = i;
			}
			lastContigID = this->entries[i].contigID;
		}
		const TotempoleEntry& lastEntry = this->entries[this->getSamples() - 1];
		contigs[lastEntry.contigID].blocksEnd = this->getSamples();
		contigs[lastEntry.contigID].maxPosition = lastEntry.maxPosition;

		// Check EOF
		reader.read(&temp_buffer[0], Constants::eof_length*sizeof(U64));
		for(U32 i = 0; i < Constants::eof_length; ++i){
			const U64* eof = reinterpret_cast<const U64*>(&temp_buffer[sizeof(U64)*i]);

			if(*eof != Constants::eof[i]){
				std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") <<  "Truncated index file!" << std::endl;
				return false;
			}
		}

		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(this->getSamples())) << " blocks..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(this->size())) << " contigs and " << Helpers::NumberThousandsSeparator(std::to_string(this->getSamples())) << " samples..." << std::endl;
		}

		U64 totalEntries = 0;
		for(U32 i = 0; i < this->getSamples(); ++i)
			totalEntries += this->entries[i].variants;

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(totalEntries)) << " variants..." << std::endl;

		// Parse
		if(!this->BuildHashTables()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TOTEMPOLE") << "Could not parse Totempole..." << std::endl;
			return false;
		}

		return true;
	}

	const U32& size(void) const{ return this->n_contigs; }
	const TotempoleEntry& operator[](const U32 p) const{ return this->entries[p]; }

	// Find overlaps function using Totempole data
	std::vector<U32> findOverlaps(const Interval& interval) const{
		std::vector<U32> ret;
		for(U32 i = this->contigs[interval.contigID].blocksStart; i < this->contigs[interval.contigID].blocksEnd; ++i){
			const TotempoleEntry& current = (*this)[i];
			if((interval.from >= current.minPosition && interval.from <= current.maxPosition) ||
					(interval.to >= current.minPosition && interval.to <= current.maxPosition) ||
					(interval.from <= current.minPosition && interval.to >= current.maxPosition))
				ret.push_back(i);

			// No need to continue searching as file is ordered

			if(current.minPosition > interval.to){
				std::cerr << "break: " << current.minPosition << ">" << interval.to << std::endl;
				break;
			}


		}

		return ret;
	}

	const U32& getLargestBlockSize(void) const{ return this->header.largest_uncompressed; }
	const U32& getBlocks(void) const{ return this->header.blocks; }
	const U64& getSamples(void) const{ return this->header.samples; }
	const contig_type& getContig(const U32 contigID) const{ return this->contigs[contigID]; }
	const header_type& getHeader(void) const{ return(this->header); }

private:
	inline bool Validate(std::ifstream& in) const{
		char MAGIC[Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH];
		in.read(MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH);

		if(strncmp(MAGIC, Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH) == 0)
			return true;
		return false;
	}

	bool BuildHashTables(void){
		if(this->size() < 1024)
			this->contigsHashTable = new hashtable(1024);
		else
			this->contigsHashTable = new hashtable(this->size() * 2);

		U32* retValue = 0;
		for(U32 i = 0; i < this->size(); ++i){
			if(this->contigsHashTable->GetItem(&this->contigs[i].name[0], &this->contigs[i].name, retValue, this->contigs[i].name.size())){
				std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Duplicated contig! Impossible!" << std::endl;
				return false;
			}
			this->contigsHashTable->SetItem(&this->contigs[i].name[0], &this->contigs[i].name, i, this->contigs[i].name.size());
		}

		if(this->getSamples() < 1024)
			this->sampleHashTable = new hashtable(1024);
		else
			this->sampleHashTable = new hashtable(this->getSamples() * 2);

		retValue = 0;
		for(U32 i = 0; i < this->getSamples(); ++i){
			if(this->sampleHashTable->GetItem(&this->samples[i][0], &this->samples[i], retValue, this->samples[i].size())){
				std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Duplicated name! Impossible!" << std::endl;
				return false;
			}
			this->sampleHashTable->SetItem(&this->samples[i][0], &this->samples[i], i, this->samples[i].size());
		}

		return true;
	}

private:
	std::string filename;
	U32 filesize;

	U32 n_contigs;
	header_type header;

	contig_type* contigs;
	std::string* samples;
	entry_type* entries;
	hashtable* contigsHashTable;
	hashtable* sampleHashTable;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOTEMPOLEREADER_H_ */
