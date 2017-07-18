#ifndef TOMAHAWK_TOTEMPOLEREADER_H_
#define TOMAHAWK_TOTEMPOLEREADER_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "MagicConstants.h"
#include "../TypeDefinitions.h"
#include "../helpers.h"
#include "TotempoleEntry.h"
#include "../algorithm/OpenHashTable.h"

namespace Tomahawk {

// Todo: temp interval
struct Interval{
	Interval(U32 contig, U32 from, U32 to) : contigID(contig), from(from), to(to){}
	~Interval(){}

	U32 contigID;
	U32 from;
	U32 to;
};

struct TotempoleContig{
typedef TotempoleContig self_type;

public:
	TotempoleContig() : minPosition(0), maxPosition(0), blocksStart(0), blocksEnd(0), bases(0){}
	~TotempoleContig(){}

	friend std::ostream& operator<<(std::ostream& stream, const TotempoleContig& entry){
		stream << entry.name << '\t' << entry.bases << '\t' << entry.minPosition << "-" << entry.maxPosition << '\t' << entry.blocksStart << "->" << entry.blocksEnd;
		return stream;
	}

	// Data is written literally to disk so needs to be reinterpreted back
	void updateFirst(char* buffer, std::ifstream& stream){
		this->bases = *reinterpret_cast<const U32*>(&buffer[0]);
		const U32* n_char = reinterpret_cast<const U32*>(&buffer[sizeof(U32)]);
		stream.read(&buffer[sizeof(U32)*2], *n_char);
		this->name = std::string(&buffer[sizeof(U32)*2], *n_char);
	}

	void updateSecond(char* buffer, std::ifstream& stream){
		this->minPosition = *reinterpret_cast<const U32*>(&buffer[0]);
		this->maxPosition = *reinterpret_cast<const U32*>(&buffer[sizeof(U32)]);
		this->blocksStart = *reinterpret_cast<const U32*>(&buffer[sizeof(U32)*2]);
		this->blocksEnd   = *reinterpret_cast<const U32*>(&buffer[sizeof(U32)*3]);
		this->bases       = *reinterpret_cast<const U32*>(&buffer[sizeof(U32)*4]);
	}

	// Updated second when read
	// contigID is implicit
	U32 minPosition;	// start position of contig
	U32 maxPosition;	// end position of contig
	U32 blocksStart; 	// start IO-seek position of blocks
	U32 blocksEnd;		// end IO-seek position of blocks

	// Updated first when read
	U32 bases;			// length of contig
	std::string name;	// contig name
};

// size of TotempoleHeader struct
#define TOTEMPOLE_headerSIZE	sizeof(float) + sizeof(U64) + sizeof(BYTE) + 3*sizeof(U32);

#pragma pack(1)
struct TotempoleHeader{
	typedef TotempoleHeader self_type;

public:
	TotempoleHeader() :
		version(0),
		samples(0),
		controller(0),
		blocks(0),
		largest_uncompressed(0),
		offset(0)
	{}

	// This ctor is used during construction
	// only possible when sample count is known
	TotempoleHeader(const U64 samples) :
		version(Constants::PROGRAM_VERSION),
		samples(samples),
		controller(0),
		blocks(0),
		largest_uncompressed(0),
		offset(0)
	{}
	~TotempoleHeader(){}

	friend std::ostream& operator<<(std::ofstream& stream, const self_type& header){
		///////////////
		// Totempole
		///////////////
		// MAGIC | version | sample count | controller byte | blocks | offset
		stream.write(Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH);
		stream.write(reinterpret_cast<const char*>(&header.version), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&header.samples), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&header.controller), sizeof(BYTE)); // Controller byte
		// At end-of-file, reopen file as in | out | binary and seek to this position and overwrite with the correct position
		stream.write(reinterpret_cast<const char*>(&header.blocks), sizeof(U32)); // Number of blocks in Tomahawk
		stream.write(reinterpret_cast<const char*>(&header.largest_uncompressed), sizeof(U32)); // Size of largest uncompressed block
		return(stream);
	}

	friend std::ostream& operator<<(std::ostream& os, const self_type& block){
		os <<
		"version: " << block.version << '\n' <<
		"samples: " << block.samples << '\n' <<
		"controller: " << std::bitset<16>(block.controller) << '\n' <<
		"blocks: " << block.blocks << '\n' <<
		"largest: " << block.largest_uncompressed << '\n' <<
		"offset: " << block.offset;

		return(os);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& block){
		stream.read(reinterpret_cast<char *>(&block.version), sizeof(float));
		stream.read(reinterpret_cast<char *>(&block.samples), sizeof(U64));
		stream.read(reinterpret_cast<char *>(&block.controller), sizeof(BYTE));
		stream.read(reinterpret_cast<char *>(&block.blocks), sizeof(U32));
		stream.read(reinterpret_cast<char *>(&block.largest_uncompressed), sizeof(U32));
		stream.read(reinterpret_cast<char *>(&block.offset), sizeof(U32));

		return(stream);
	}

public:
	float version;		// version used to write header
	U64 samples;		// number of samples
	BYTE controller;	// controller block
	U32 blocks;			// number of blocks in Tomahawk
	U32 largest_uncompressed;	// largest block-size in bytes
	U32 offset;			// IO disk offset for start of data
};

class TotempoleReader {
	typedef TotempoleHeader headertype;
	typedef TotempoleContig contig_type;
	typedef TotempoleEntry  entry_type;

	typedef Tomahawk::Hash::HashTable<std::string, U32> hashtable;

public:
	TotempoleReader() : filesize(0), contigs(nullptr), samples(nullptr), entries(nullptr), contigsHashTable(nullptr), sampleHashTable(nullptr){}
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

		// Get number of contigs
		reader.read(reinterpret_cast<char *>(&this->n_contigs), sizeof(U32));

		char temp_buffer[65536];
		this->contigs = new contig_type[this->size()];
		for(U32 i = 0; i < this->size(); ++i){
			reader.read(&temp_buffer[0], sizeof(U32)*2);
			this->contigs[i].updateFirst(temp_buffer, reader);
			//std::cerr << contigs[i] << std::endl;
		}


		this->samples = new std::string[this->header.samples];
		for(U32 i = 0; i < this->header.samples; ++i){
			reader.read(&temp_buffer[0], sizeof(U32));
			const U32 length = *reinterpret_cast<const U32*>(&temp_buffer[0]);
			reader.read(&temp_buffer[sizeof(U32)], length);
			this->samples[i] = std::string(&temp_buffer[sizeof(U32)], length);
			//std::cerr << i << '\t' << samples[i] << std::endl;
		}

		if(reader.tellg() != this->header.offset){
			std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Corrupt file" << std::endl;
			std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << reader.tellg() << '/' << this->header.offset << std::endl;
			return false;
		}

		// Populate Totempole entries
		this->entries = new TotempoleEntry[this->header.blocks];
		for(U32 i = 0; i < this->header.blocks; ++i){
			reader >> this->entries[i];
			//std::cerr << this->entries[i] << std::endl;
		}

		// Find boundaries for Totempole blocks
		// Master index of indices
		U32 lastContigID = this->entries[0].contigID;
		contigs[lastContigID].minPosition = this->entries[0].minPosition;
		contigs[lastContigID].blocksStart = 0;
		for(U32 i = 1; i < this->header.blocks; ++i){
			if(lastContigID != this->entries[i].contigID){
				contigs[lastContigID].maxPosition = this->entries[i-1].maxPosition;
				contigs[lastContigID].blocksEnd = i;
				contigs[this->entries[i].contigID].minPosition = this->entries[i].minPosition;
				contigs[this->entries[i].contigID].blocksStart = i;
			}
			lastContigID = this->entries[i].contigID;
		}
		const TotempoleEntry& lastEntry = this->entries[this->header.blocks - 1];
		contigs[lastEntry.contigID].blocksEnd = this->header.blocks;
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
			std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(this->header.blocks)) << " blocks..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(this->size())) << " contigs and " << Helpers::NumberThousandsSeparator(std::to_string(this->header.samples)) << " samples..." << std::endl;
		}

		U64 totalEntries = 0;
		for(U32 i = 0; i < this->header.blocks; ++i)
			totalEntries += this->entries[i].variants;

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(totalEntries)) << " variants..." << std::endl;

		// Parse
		if(!this->Parse()){
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
	const TotempoleContig& getContig(const U32 contigID) const{ return this->contigs[contigID]; }
	const TotempoleHeader& getHeader(void) const{ return(this->header); }

private:
	inline bool Validate(std::ifstream& in) const{
		char MAGIC[Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH];
		in.read(MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH);

		if(strncmp(MAGIC, Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH) == 0)
			return true;
		return false;
	}

	bool Parse(){
		if(!this->BuildHashTables())
			return false;

		return true;
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

		if(this->header.samples < 1024)
			this->sampleHashTable = new hashtable(1024);
		else
			this->sampleHashTable = new hashtable(this->header.samples * 2);

		retValue = 0;
		for(U32 i = 0; i < this->header.samples; ++i){
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
	headertype header;

	contig_type* contigs;
	std::string* samples;
	entry_type* entries;
	hashtable* contigsHashTable;
	hashtable* sampleHashTable;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOTEMPOLEREADER_H_ */
