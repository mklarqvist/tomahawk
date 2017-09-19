#include "TotempoleOutputReader.h"

namespace Tomahawk {
namespace Totempole {

TotempoleOutputReader::TotempoleOutputReader() : n_entries(0), ERROR_STATE(TOI_OK), entries(nullptr){}
TotempoleOutputReader::~TotempoleOutputReader(){}

bool TotempoleOutputReader::Open(const std::string& input){
	this->stream.open(input, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream.good()){
		//std::cerr << "failed does not exist" << std::endl;
		this->ERROR_STATE = TOI_NO_EXIST;
		return false;
	}

	U32 filesize = this->stream.tellg();
	this->stream.seekg(0);

	this->stream >> this->header;
	if(!this->header.validate(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC)){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "Incorrect header" << std::endl;
		//std::cerr << this->header.
		this->ERROR_STATE = TOI_CORRUPTED;
		exit(1);
	}

	this->buffer.resize(filesize);
	const U32 readUntil = filesize - (U32)this->stream.tellg();

	if(readUntil % sizeof(entry_type) != 0){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "Corrupted data!" << std::endl;
		this->ERROR_STATE = TOI_CORRUPTED;
		exit(1);
	}

	this->stream.read(this->buffer.data, readUntil);
	this->entries = reinterpret_cast<const entry_type*>(this->buffer.data);
	this->n_entries = readUntil / sizeof(entry_type);

	//for(U32 i = 0; i < this->n_entries; ++i)
	//	std::cerr << this->entries[i] << std::endl;

	this->ERROR_STATE = TOI_OK;
	return true;
}

bool TotempoleOutputReader::findOverlap(const S32 contigID){
	if(this->ERROR_STATE != TOI_OK){
		std::cerr << "no toi" << std::endl;
		return false;
	}

	if(!this->header.controller.sorted){
		std::cerr << "controller not sorted" << std::endl;
		return false;
	}

	std::cerr << "searching for: " << contigID << std::endl;
	for(U32 i = 0; i < this->size(); ++i){

	}

	return true;
}


} /* namespace Totempole */
} /* namespace Tomahawk */
