#include "TotempoleReader.h"

namespace Tomahawk {

TotempoleReader::TotempoleReader() :
		filesize(0),
		n_contigs(0),
		contigs(nullptr),
		samples(nullptr),
		entries(nullptr),
		contigsHashTable(nullptr),
		sampleHashTable(nullptr)
{}

TotempoleReader::~TotempoleReader(){
	delete [] this->contigs;
	delete [] this->entries;
	delete [] this->samples;

	delete this->contigsHashTable;
	delete this->sampleHashTable;
}

bool TotempoleReader::Validate(std::ifstream& in) const{
	char MAGIC[Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH];
	in.read(MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH);

	if(strncmp(MAGIC, Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH) == 0)
		return true;
	return false;
}

bool TotempoleReader::Open(const std::string filename){
	if(filename.size() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "IO") << "No input filename..." << std::endl;
		return false;
	}

	this->filename = filename;

	this->stream.open(this->filename, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "IO") << "Could not open: " << this->filename << "..." << std::endl;
		return false;
	}

	this->filesize = this->stream.tellg();
	this->stream.seekg(0);

	if(this->filesize <= 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "IO") << "File size is 0..." << std::endl;
		return false;
	}


	if(this->filesize < Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TOTEMPOLE") << "Failed MAGIC..." << std::endl;
		return false;
	}

	// Reader header and validate
	if(!this->Validate(this->stream)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TOTEMPOLE") << "Could not validate Totempole header..." << std::endl;
		return false;
	}
	// Load header data
	this->stream >> this->header;
#if DEBUG_MODE == 1
	std::cerr << this->header << std::endl;
#endif


	// Get number of contigs
	this->stream.read(reinterpret_cast<char *>(&this->n_contigs), sizeof(U32));
#if DEBUG_MODE == 1
	std::cerr << this->n_contigs << std::endl;
#endif
	this->contigs = new contig_type[this->size()];
	for(U32 i = 0; i < this->size(); ++i){
		contig_base_type* contig_base = reinterpret_cast<contig_base_type*>(&this->contigs[i]);
		this->stream >> *contig_base;
#if DEBUG_MODE == 1
		std::cerr << *contig_base << std::endl;
#endif
	}

	char temp_buffer[65536];
	this->samples = new std::string[this->getSamples()];
	for(U32 i = 0; i < this->getSamples(); ++i){
		this->stream.read(&temp_buffer[0], sizeof(U32));
		const U32 length = *reinterpret_cast<const U32*>(&temp_buffer[0]);
		this->stream.read(&temp_buffer[sizeof(U32)], length);
		this->samples[i] = std::string(&temp_buffer[sizeof(U32)], length);
#if DEBUG_MODE == 1
		std::cerr << i << '\t' << samples[i] << std::endl;
#endif
	}

	if(this->stream.tellg() != this->header.offset){
		std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Corrupt file" << std::endl;
		std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << this->stream.tellg() << '/' << this->header.offset << std::endl;
		return false;
	}

#if DEBUG_MODE == 1
	std::cerr << this->n_contigs << '\t' << this->header.blocks << '\t' << this->header.samples << std::endl;
#endif

	// Populate Totempole entries
	this->entries = new entry_type[this->getBlocks()];
	for(U32 i = 0; i < this->getBlocks(); ++i){
		this->stream >> this->entries[i];
#if DEBUG_MODE == 1
		std::cerr << i << '\t' << this->header.blocks << '\t' << this->entries[i] << std::endl;
#endif
	}

	this->BuildUpdateContigs();

#if DEBUG_MODE == 1
	for(U32 i = 0; i < this->size(); ++i)
		std::cerr << this->contigs[i] << std::endl;
#endif

	if(!this->ValidateEOF(this->stream))
		return false;

	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(this->getSamples())) << " blocks..." << std::endl;
		std::cerr << Helpers::timestamp("LOG", "TOTEMPOLE") << "Found: " << Helpers::NumberThousandsSeparator(std::to_string(this->size())) << " contigs and " << Helpers::NumberThousandsSeparator(std::to_string(this->getSamples())) << " samples..." << std::endl;
	}

	U64 totalEntries = 0;
	for(U32 i = 0; i < this->getBlocks(); ++i)
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

bool TotempoleReader::ValidateEOF(std::ifstream& in){
	char temp_buffer[Constants::eof_length*sizeof(U64)];
	in.read(&temp_buffer[0], Constants::eof_length*sizeof(U64));
	for(U32 i = 0; i < Constants::eof_length; ++i){
		const U64* eof = reinterpret_cast<const U64*>(&temp_buffer[sizeof(U64)*i]);

		if(*eof != Constants::eof[i]){
			std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") <<  "Truncated index file!" << std::endl;
			return false;
		}
	}
	return true;
}

void TotempoleReader::BuildUpdateContigs(void){
	// Find boundaries for Totempole blocks
	// Master index of indices
	// Update contig data
	U32 lastContigID = this->entries[0].contigID;
	this->contigs[lastContigID].minPosition = this->entries[0].minPosition;
	this->contigs[lastContigID].blocksStart = 0;
	for(U32 i = 1; i < this->getBlocks(); ++i){
		if(lastContigID != this->entries[i].contigID){
			this->contigs[lastContigID].maxPosition = this->entries[i-1].maxPosition;
			this->contigs[lastContigID].blocksEnd = i;
			this->contigs[this->entries[i].contigID].minPosition = this->entries[i].minPosition;
			this->contigs[this->entries[i].contigID].blocksStart = i;
		}
		lastContigID = this->entries[i].contigID;
	}
	const TotempoleEntry& lastEntry = this->entries[this->getBlocks() - 1];
	this->contigs[lastEntry.contigID].blocksEnd = this->getBlocks();
	this->contigs[lastEntry.contigID].maxPosition = lastEntry.maxPosition;
}

bool TotempoleReader::BuildHashTables(void){
	if(this->size() < 1024)
		this->contigsHashTable = new hash_table(1024);
	else
		this->contigsHashTable = new hash_table(this->size() * 2);

	U32* retValue = 0;
	for(U32 i = 0; i < this->size(); ++i){
		if(this->contigsHashTable->GetItem(&this->contigs[i].name[0], &this->contigs[i].name, retValue, this->contigs[i].name.size())){
			std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Duplicated contig! Impossible!" << std::endl;
			return false;
		}
		this->contigsHashTable->SetItem(&this->contigs[i].name[0], &this->contigs[i].name, i, this->contigs[i].name.size());
	}

	if(this->getSamples() < 1024)
		this->sampleHashTable = new hash_table(1024);
	else
		this->sampleHashTable = new hash_table(this->getSamples() * 2);

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

// Find overlaps function using Totempole data
std::vector<U32> TotempoleReader::findOverlaps(const Interval& interval) const{
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

} /* namespace Tomahawk */
