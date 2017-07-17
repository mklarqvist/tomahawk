#ifndef TOMAHAWK_TOTEMPOLEREADER_H_
#define TOMAHAWK_TOTEMPOLEREADER_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "MagicConstants.h"
#include "../TypeDefinitions.h"
#include "../helpers.h"
#include "TotempoleEntry.h"
#include "../algorithm/OpenHashTable.h"

namespace Tomahawk {

struct Interval{
	Interval(U32 contig, U32 from, U32 to) : contigID(contig), from(from), to(to){}
	~Interval(){}

	U32 contigID;
	U32 from;
	U32 to;
};

struct TotempoleContig{
	TotempoleContig() : minPosition(0), maxPosition(0), blocksStart(0), blocksEnd(0), bases(0){}
	~TotempoleContig(){}

	friend std::ostream& operator<<(std::ostream& stream, const TotempoleContig& entry){
		stream << entry.name << '\t' << entry.bases << '\t' << entry.minPosition << "-" << entry.maxPosition << '\t' << entry.blocksStart << "->" << entry.blocksEnd;
		return stream;
	}

	bool operator>>(std::istream& stream); // Overload when reading from cin or ifstream
	//friend std::ofstream& operator<<();

	// contigID is implicit
	U32 minPosition;	// start position of contig
	U32 maxPosition;	// end position of contig
	U32 blocksStart; 	// start IO-seek position of blocks
	U32 blocksEnd;		// end IO-seek position of blocks
	U32 bases;			// length of contig
	std::string name;	// contig name
};

namespace Constants{
const U32 TOTEMPOLE_HEADER_SIZE = sizeof(float) + sizeof(U64) + sizeof(BYTE) + 3*sizeof(U32);
}

#pragma pack(1)
struct TotempoleHeader{
public:
	TotempoleHeader() :
		version(Constants::PROGRAM_VERSION),
		samples(0),
		controller(0),
		blocks(0),
		largest_uncompressed(0),
		offset(0)
	{}
	TotempoleHeader(const U64 samples) :
		version(Constants::PROGRAM_VERSION),
		samples(samples),
		controller(0),
		blocks(0),
		largest_uncompressed(0),
		offset(0)
	{}
	~TotempoleHeader(){}

	friend std::ostream& operator<<(std::ofstream& stream, const TotempoleHeader& header){
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

public:
	float version;		// version used to write header
	U64 samples;		// number of samples
	BYTE controller;	// controller block
	U32 blocks;			// number of blocks in Tomahawk
	U32 largest_uncompressed;	// largest block-size in bytes
	U32 offset;			// IO disk offset for start of data
};

class TotempoleReader {
	typedef Tomahawk::Hash::HashTable<std::string, U32> hashtable;

public:
	TotempoleReader() : filesize_(0), data_(nullptr), header_(nullptr), entries_(nullptr), contigsHashTable_(nullptr), sampleHashTable_(nullptr){}
	~TotempoleReader(){
		delete [] this->data_;
		delete this->contigsHashTable_;
		delete this->sampleHashTable_;
	}

	bool Open(const std::string filename){
		if(filename.size() == 0){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "IO") << "No input filename..." << std::endl;
			return false;
		}

		this->filename_ = filename;

		std::ifstream reader(this->filename_, std::ios::in | std::ios::binary | std::ios::ate);
		if(!reader.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "IO") << "Could not open: " << this->filename_ << "..." << std::endl;
			return false;
		}

		this->filesize_ = reader.tellg();
		reader.seekg(0);

		if(this->filesize_ <= 0){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "IO") << "File size is 0..." << std::endl;
			return false;
		}

		this->data_ = new char[this->filesize_];
		reader.read(this->data_, this->filesize_);
		reader.close();

		// Parse
		if(!this->Parse()){
			std::cerr << "Could not parse" << std::endl;
			return false;
		}

		return true;
	}

	U32 size(void) const{ return this->contigs_.size(); }
	const TotempoleEntry& operator[](const U32 p) const{ return this->entries_[p]; }

	// Find overlaps function using Totempole data
	std::vector<U32> findOverlaps(const Interval& interval) const{
		std::vector<U32> ret;
		for(U32 i = this->contigs_[interval.contigID].blocksStart; i < this->contigs_[interval.contigID].blocksEnd; ++i){
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

	const U32& getLargestBlockSize(void) const{ return this->header_->largest_uncompressed; }
	const U32& getBlocks(void) const{ return this->header_->blocks; }
	const U64& getSamples(void) const{ return this->header_->samples; }
	const TotempoleContig& getContig(const U32 contigID) const{ return this->contigs_[contigID]; }
	const TotempoleHeader& getHeader(void) const{ return(*this->header_); }

private:
	inline bool Validate() const{
		if(strncmp(this->data_, Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH) == 0)
			return true;
		return false;
	}

	bool Parse(){
		if(this->filesize_ < Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH){
			std::cerr << "file too small" << std::endl;
			return false;
		}

		if(!this->Validate()){
			std::cerr << "could not vlaidate" << std::endl;
			return false;
		}

		U32 buffer_pointer = Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH;
		this->header_ = reinterpret_cast<TotempoleHeader*>(&this->data_[buffer_pointer]);
		buffer_pointer += Constants::TOTEMPOLE_HEADER_SIZE;

		std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(this->header_->blocks)) << " blocks..." << std::endl;

		U32* contigs = reinterpret_cast<U32*>(&this->data_[buffer_pointer]);
		buffer_pointer += sizeof(U32);

		for(U32 i = 0; i < *contigs; ++i){
			TotempoleContig contig;
			U32* bases = reinterpret_cast<U32*>(&this->data_[buffer_pointer]);
			buffer_pointer += sizeof(U32);
			contig.bases = *bases;
			U32* char_length = reinterpret_cast<U32*>(&this->data_[buffer_pointer]);
			buffer_pointer += sizeof(U32);
			contig.name = std::string(&this->data_[buffer_pointer], *char_length);
			this->contigs_.push_back(contig);
			buffer_pointer += *char_length;
		}

		for(U32 i = 0; i < this->header_->samples; ++i){
			U32* length = reinterpret_cast<U32*>(&this->data_[buffer_pointer]);
			buffer_pointer += sizeof(U32);
			this->sampleNames_.push_back(std::string(&this->data_[buffer_pointer], *length));
			buffer_pointer += *length;
		}

		std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(this->contigs_.size())) << " contigs and " << Helpers::NumberThousandsSeparator(std::to_string(this->sampleNames_.size())) << " samples..." << std::endl;

		if(buffer_pointer != this->header_->offset){
			std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Corrupt file" << std::endl;
			std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << buffer_pointer << '/' << this->header_->offset << std::endl;
			return false;
		}

		this->entries_ = reinterpret_cast<TotempoleEntry*>(&this->data_[buffer_pointer]);
		buffer_pointer += this->header_->blocks * Constants::TOTEMPOLE_ENTRY_SIZE;

		// Parse
		U32 lastContigID = this->entries_[0].contigID;
		this->contigs_[lastContigID].minPosition = this->entries_[0].minPosition;
		this->contigs_[lastContigID].blocksStart = 0;
		for(U32 i = 1; i < this->header_->blocks; ++i){
			if(lastContigID != this->entries_[i].contigID){
				this->contigs_[lastContigID].maxPosition = this->entries_[i-1].maxPosition;
				this->contigs_[lastContigID].blocksEnd = i;
				this->contigs_[this->entries_[i].contigID].minPosition = this->entries_[i].minPosition;
				this->contigs_[this->entries_[i].contigID].blocksStart = i;
			}
			lastContigID = this->entries_[i].contigID;
		}
		const TotempoleEntry& lastEntry = this->entries_[this->header_->blocks - 1];
		this->contigs_[lastEntry.contigID].blocksEnd = this->header_->blocks;
		this->contigs_[lastEntry.contigID].maxPosition = lastEntry.maxPosition;

		// temp
		/*
		for(U32 i = 0; i < this->contigs_.size(); ++i){
			std::cerr << this->contigs_[i] << std::endl;
			//U32 width = this->contigs_[i].blocksEnd - this->contigs_[i].blocksStart;
			//std::cerr << "Width: " << width << std::endl;

			for(U32 j = this->contigs_[i].blocksStart; j < this->contigs_[i].blocksEnd; ++j){
				std::cerr << '\t' << (*this)[j].byte_offset << '\t' << (*this)[j].minPosition << "->" << (*this)[j].maxPosition << '\t' << (*this)[j].variants << std::endl;
			}

		}
		//
		 */

		U64 totalEntries = 0;
		for(U32 i = 0; i < this->header_->blocks; ++i)
			totalEntries += (*this)[i].variants;

		std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(totalEntries)) << " variants..." << std::endl;

		for(U32 i = 0; i < Constants::eof_length; ++i){
			U64* eof = reinterpret_cast<U64*>(&this->data_[buffer_pointer]);
			buffer_pointer += sizeof(U64);

			if(*eof != Constants::eof[i]){
				std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") <<  "Truncated index file!" << std::endl;
				return false;
			}
		}

		if(!this->BuildHashTables())
			return false;

		return true;
	}

	bool BuildHashTables(void){
		if(this->contigs_.size() < 1024)
			this->contigsHashTable_ = new hashtable(1024);
		else
			this->contigsHashTable_ = new hashtable(this->contigs_.size() * 2);

		U32* retValue = 0;
		for(U32 i = 0; i < this->contigs_.size(); ++i){
			if(this->contigsHashTable_->GetItem(&this->contigs_[i].name[0], &this->contigs_[i].name, retValue, this->contigs_[i].name.size())){
				std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Duplicated contig! Impossible!" << std::endl;
				return false;
			}
			this->contigsHashTable_->SetItem(&this->contigs_[i].name[0], &this->contigs_[i].name, i, this->contigs_[i].name.size());
		}

		if(this->sampleNames_.size() < 1024)
			this->sampleHashTable_ = new hashtable(1024);
		else
			this->sampleHashTable_ = new hashtable(this->sampleNames_.size() * 2);

		retValue = 0;
		for(U32 i = 0; i < this->sampleNames_.size(); ++i){
			if(this->sampleHashTable_->GetItem(&this->sampleNames_[i][0], &this->sampleNames_[i], retValue, this->sampleNames_[i].size())){
				std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Duplicated name! Impossible!" << std::endl;
				return false;
			}
			this->sampleHashTable_->SetItem(&this->sampleNames_[i][0], &this->sampleNames_[i], i, this->sampleNames_[i].size());
		}

		return true;
	}

private:
	std::string filename_;
	U32 filesize_;
	std::vector<TotempoleContig> contigs_;
	std::vector<std::string> sampleNames_;
	char* data_;
	TotempoleHeader* header_;
	TotempoleEntry* entries_;
	hashtable* contigsHashTable_;
	hashtable* sampleHashTable_;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOTEMPOLEREADER_H_ */
