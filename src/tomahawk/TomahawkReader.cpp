#include "TomahawkReader.h"

namespace Tomahawk {

// Remember to resize buffers to header.getLargestBlockSize()+64 after header is loaded
TomahawkReader::TomahawkReader() :
	samples(0),
	version(0),
	filesize(0),
	bit_width(0),
	dropGenotypes(false),
	showHeader(true),
	currentBlockID(0),
	writer(nullptr),
	group_htable(nullptr)
{}

TomahawkReader::~TomahawkReader(){
	this->buffer.deleteAll();
	this->data.deleteAll();
	this->outputBuffer.deleteAll();
	delete writer;
	delete group_htable;
}

bool TomahawkReader::Open(const std::string input){
	const std::string index = input + '.' + Tomahawk::Constants::OUTPUT_INDEX_SUFFIX;

	// Parse Totempole
	if(!this->totempole.Open(index)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TOTEMPOLE") << "Failed build!" << std::endl;
		return false;
	}

	// Resize buffers to accomodate the largest possible block
	// without ever resizing
	// this is for performance reasons
	this->buffer.resize(this->totempole.getLargestBlockSize() + 64);
	this->data.resize(this->totempole.getLargestBlockSize() + 64);
	this->outputBuffer.resize(this->totempole.getLargestBlockSize() + 64);

	if(input.size() == 0){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "No input filename..." << std::endl;
		return false;
	}

	this->stream.open(input, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Could not open " << input << "..." << std::endl;
		return false;
	}
	this->filesize = this->stream.tellg();
	this->stream.seekg(0);

	// Validate MAGIC and header
	if(!this->ValidateHeader(this->stream)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to validate header..." << std::endl;
		return false;
	}

	return(this->Validate());
}


inline bool TomahawkReader::ValidateHeader(std::ifstream& in) const{
	char MAGIC[Constants::WRITE_HEADER_MAGIC_LENGTH];
	in.read(MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH);

	if(strncmp(&MAGIC[0], &Constants::WRITE_HEADER_MAGIC[0], Constants::WRITE_HEADER_MAGIC_LENGTH) == 0)
		return true;

	return false;
}

bool TomahawkReader::loadGroups(const std::string& file){
	if(this->Occ.size() != 0){
		std::cerr << "groups already loaded" << std::endl;
		return false;
	}

	if(file.size() == 0){
		std::cerr << "No file set" << std::endl;
		return false;
	}

	// Open stream
	std::ifstream stream(file);
	if(!stream.good()){
		std::cerr << "bad file" << std::endl;
		return false;
	}

	this->group_htable = new hash_table(10000);


	std::string line;
	std::vector< std::vector< S32 > > groupings(this->samples);

	S32* sampleID = nullptr;
	S32* groupID_lookup = nullptr;
	S32 groupID = 0;
	U32 n_lines = 0;

	while(getline(stream, line)){
		// Empty lines
		if(line.size() == 0)
			break;

		// Assert correct format
		// Count tabs until out-of-range
		size_t prev_pos = 0;
		size_t pos = line.find('\t', prev_pos + 1);

		if(pos == std::string::npos){
			std::cerr << "illegal: has no data for sample (" << n_lines << ")" << std::endl;
			return(false);
		}

		// Make sure the sample name exists in this file
		const std::string sampleName = std::string(&line[prev_pos], pos - prev_pos);
		if(!this->totempole.sampleHashTable->GetItem(&sampleName[0], &sampleName, sampleID, sampleName.length())){
			std::cerr << "sample does not exist" << std::endl;
			return false;
		}

		// Cycle over lines
		prev_pos = pos + 1;
		while(prev_pos != std::string::npos){
			pos = line.find('\t', prev_pos + 1);
			if(pos == std::string::npos){
				pos = line.size();
			}

			const std::string group = std::string(&line[prev_pos], pos - prev_pos);
			if(!this->group_htable->GetItem(&group[0], &group, groupID_lookup, group.length())){
				this->group_htable->SetItem(&group[0], &group, groupID, group.length());
				this->groups.push_back(GroupPair(group));

				groupings[*sampleID].push_back(groupID);
				++groupID;
			} else {
				groupings[*sampleID].push_back(*groupID_lookup);
				++this->groups[*groupID_lookup];
			}

			if(pos == line.size())
				break;

			prev_pos = pos + 1;
		}

		++n_lines;
	}

	if(groupID == 0){
		std::cerr << "no data loaded" << std::endl;
		return false;
	}

	this->Occ = occ_matrix(this->samples + 1, std::vector< U64 >(groupID, 0));
	occ_vector* prev = &Occ[0];

	for(U32 i = 0; i < this->samples; ++i){
		// Propagate previous vector counts
		for(U32 j = 0; j < groupID; ++j)
			this->Occ[i + 1][j] = prev->at(j);

		// Update
		for(U32 j = 0; j < groupings[i].size(); ++j)
			this->Occ[i + 1][groupings[i][j]] = prev->at(groupings[i][j]) + 1;

		prev = &this->Occ[i + 1];
	}

	for(U32 i = 0; i < this->groups.size(); ++i)
		std::cerr << i << '\t' << this->groups[i].name << '\t' << this->groups[i].n_entries << std::endl;

	// Temp
	// Dump
	/*
	for(U32 i = 0; i < this->Occ.size(); ++i){
		for(U32 j = 0; j < this->Occ[0].size(); ++j){
			std::cout << this->Occ[i][j] << '\t';
		}
		std::cout << std::endl;
	}
	*/

	return true;
}

bool TomahawkReader::getBlocks(void){
	U64 buffer_size = 0;
	for(U32 i = 0; i < this->totempole.header.blocks; ++i){
		buffer_size += 0;
	}

	if(buffer_size == 0){
		std::cerr << "impsosible" << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << this->totempole.getHeader().blocks << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

	this->data.resize(buffer_size + 1000);

	for(U32 i = 0; i < this->totempole.getBlocks(); ++i)
		if(!this->getBlock(i))
			return false;

	return true;
}

bool TomahawkReader::getBlocks(std::vector<U32>& blocks){
	U64 buffer_size = 0;
	for(U32 i = 0; i < blocks.size(); ++i){
		//std::cerr << i << '/' << blocks.size() << '\t' << this->totempole[i].l_uncompressed << std::endl;
		buffer_size += 0;
	}

	if(buffer_size == 0)
		return false;

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << blocks.size() << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

	this->data.resize(buffer_size + 1000);

	for(U32 i = 0; i < blocks.size(); ++i)
		if(!this->getBlock(blocks[i]))
			return false;

	return true;
}

bool TomahawkReader::getBlocks(std::vector< std::pair<U32, U32> >& blocks){
	U64 buffer_size = 0;
	for(U32 i = 0; i < blocks.size(); ++i){
		for(U32 j = blocks[i].first; j < blocks[i].second; ++j){
			buffer_size += 0;
		}
	}

	if(buffer_size == 0){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Uncompressed size is 0..." << std::endl;
		return false;
	}

	U64 totalBlocks = 0;
	for(U32 i = 0; i < blocks.size(); ++i)
		totalBlocks += blocks[i].second - blocks[i].first;

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflating " << totalBlocks << " blocks into " << Helpers::ToPrettyString(buffer_size/1000) << " kb..." << std::endl;

	this->data.resize(buffer_size + 1000);

	for(U32 i = 0; i < blocks.size(); ++i){
		for(U32 j = blocks[i].first; j < blocks[i].second; ++j){
			if(!this->getBlock(j)){
				std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Could not get block: " << j << '/' << blocks[i].second << std::endl;
				return false;
			}
		}
	}
	return true;
}


bool TomahawkReader::getBlock(const U32 blockID){
	if(!this->stream.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Bad file stream for block " << blockID << "..." << std::endl;
		return false;
	}

	this->stream.seekg(this->totempole[blockID].byte_offset);
	if(!this->stream.good()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed search..." << std::endl;
		return false;
	}

	const U32 readLength = this->totempole[blockID].byte_offset_end - this->totempole[blockID].byte_offset;

	if(readLength > this->buffer.capacity()){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Overflowing capacity: " << readLength << '/' << this->buffer.capacity() << std::endl;
		exit(1);
	}

	this->blockDataOffsets.push_back(DataOffsetPair(&this->data.data[this->data.pointer], this->totempole[blockID]));
	if(!this->stream.read(&this->buffer.data[0], readLength)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed read: " << this->stream.good() << '\t' << this->stream.fail() << '/' << this->stream.eof() << std::endl;
		//std::cerr << this->stream_.gcount() << '/' << readLength << std::endl;
		return false;
	}
	this->buffer.pointer = readLength;

	if(!this->tgzf_controller.Inflate(this->buffer, this->data)){
		std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to inflate data..." << std::endl;
		return false;
	}

	return true;
}

bool TomahawkReader::Validate(void){
	char temp_buffer[sizeof(float)+sizeof(U64)];
	this->stream.read(&temp_buffer[0], sizeof(float)+sizeof(U64));

	const float* version = reinterpret_cast<const float*>(&temp_buffer[0]);
	const U64* samples   = reinterpret_cast<const U64*>(&temp_buffer[sizeof(float)]);

	this->version = *version;
	this->samples = *samples;

	if(this->version != this->totempole.getHeader().version){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "File discordance: versions do not match..." << std::endl;
		return false;
	}

	if(this->samples != this->totempole.getHeader().samples){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "File discordance:number of samples do not match" << std::endl;
		return false;
	}

	// Determine what bit-width functions to use
	this->DetermineBitWidth();

	return true;
}

void TomahawkReader::DetermineBitWidth(void){
	if(this->samples <= Constants::UPPER_LIMIT_SAMPLES_8B){
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_8B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 8-bit width..." << std::endl;
		}
		this->bit_width = sizeof(BYTE);
	} else if(this->samples <= Constants::UPPER_LIMIT_SAMPLES_16B){
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_16B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 16-bit width..." << std::endl;
		}
		this->bit_width = sizeof(U16);
	} else if(this->samples <= Constants::UPPER_LIMIT_SAMPLES_32B){
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_32B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 32-bit width..." << std::endl;
		}
		this->bit_width = sizeof(U32);
	} else {
		if(!SILENT){
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_8B  << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_16B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " > " << Constants::UPPER_LIMIT_SAMPLES_32B << "... Skip" << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Samples: " << this->samples << " < " << Constants::UPPER_LIMIT_SAMPLES_64B << "..." << std::endl;
			std::cerr << Helpers::timestamp("LOG", "RLE") << "Using 64-bit width..." << std::endl;
		}
		this->bit_width = sizeof(U64);
	}
}

bool TomahawkReader::outputBlocks(std::vector<U32>& blocks){
	typedef bool (Tomahawk::TomahawkReader::*blockFunction)(const U32& blockID);
	blockFunction func__ = nullptr;

	switch(this->bit_width){
	case 1: func__ = &TomahawkReader::__outputBlock<BYTE>; break;
	case 2: func__ = &TomahawkReader::__outputBlock<U16>;  break;
	case 4: func__ = &TomahawkReader::__outputBlock<U32>;  break;
	case 8: func__ = &TomahawkReader::__outputBlock<U64>;  break;
	default: exit(1); break;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","TGZF") << "Inflating " << blocks.size() << " blocks..." << std::endl;

	// Output header
	if(this->showHeader)
		std::cout << this->totempole.literals + "\n##tomahawk_viewCommand=" + Helpers::program_string() << std::endl;

	for(U32 i = 0; i < blocks.size(); ++i){
		this->getBlock(blocks[i]);
		(*this.*func__)(blocks[i]);
	}

	return true;
}

bool TomahawkReader::outputBlocks(){
	typedef bool (Tomahawk::TomahawkReader::*blockFunction)(const U32& blockID);
	blockFunction func__ = nullptr;

	switch(this->bit_width){
	case 1: func__ = &TomahawkReader::__outputBlock<BYTE>; break;
	case 2: func__ = &TomahawkReader::__outputBlock<U16>;  break;
	case 4: func__ = &TomahawkReader::__outputBlock<U32>;  break;
	case 8: func__ = &TomahawkReader::__outputBlock<U64>;  break;
	default:
		std::cerr << Helpers::timestamp("ERROR") << "Word sizing could not be determined!" << std::endl;
		exit(1);
		break;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "TGZF") << "Inflating " << this->totempole.getBlocks() << " blocks..." << std::endl;

	// Output header
	if(this->showHeader){
		std::cout << this->totempole.literals << std::endl;
		std::cout << "##INFO=<ID=HWE_P,Number=1,Type=Float,Description=\"Hardy-Weinberg P-value (Fisher's exact test)\">" << std::endl;
		std::cout << "##INFO=<ID=MGF,Number=1,Type=Float,Description=\"Minor genotype frequency\">" << std::endl;
		std::cout << "##tomahawk_viewCommand=" + Helpers::program_string() << std::endl;
		std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
		for(U32 i = 0; i < this->totempole.header.samples - 1; ++i)
			std::cout << this->totempole.samples[i] << '\t';
		std::cout << this->totempole.samples[this->totempole.header.samples - 1] << std::endl;
	}

	for(U32 i = 0; i < this->totempole.getBlocks(); ++i){
		this->getBlock(i);
		(*this.*func__)(i);
	}

	return true;
}

bool TomahawkReader::nextBlock(const bool clear){
	if(!this->stream.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Stream bad " << this->currentBlockID << std::endl;
		return false;
	}

	this->buffer.reset();
	this->data.reset();

	const U32 readLength = this->totempole[this->currentBlockID].byte_offset_end - this->totempole[this->currentBlockID].byte_offset;

	// Read from start to start + byte-width
	if(!this->stream.read(&this->buffer.data[0], readLength)){
		std::cerr << Helpers::timestamp("ERROR", "TWK") << "Failed read: " << this->stream.good() << '\t' << this->stream.fail() << '/' << this->stream.eof() << std::endl;
		std::cerr << this->stream.gcount() << '/' << readLength << std::endl;
		return false;
	}

	// Set buffer width to data loaded size
	this->buffer.pointer = readLength;

	// Inflate TGZF block
	if(!this->tgzf_controller.Inflate(this->buffer, this->data)){
		std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to inflate DATA..." << std::endl;
		return false;
	}

	// Update blockID
	++this->currentBlockID;

	return true;
}

} /* namespace Tomahawk */
