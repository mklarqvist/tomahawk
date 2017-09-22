#include "TotempoleOutputReader.h"

namespace Tomahawk {
namespace Totempole {

TotempoleOutputReader::TotempoleOutputReader() : ERROR_STATE(TOI_INIT), n_entries(0), entries(nullptr), index(nullptr){}
TotempoleOutputReader::~TotempoleOutputReader(){
	delete this->index;
}

bool TotempoleOutputReader::Open(const std::string& input, const contig_type* contigs){
	this->stream.open(input, std::ios::in | std::ios::binary | std::ios::ate);
	if(!this->stream.good()){
		//std::cerr << "failed does not exist" << std::endl;
		this->ERROR_STATE = TOI_NO_EXIST;
		return false;
	}

	U64 filesize = this->stream.tellg();
	this->stream.seekg(0);

	this->stream >> this->header;
	if(!this->header.validate(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC)){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "Incorrect header!" << std::endl;
		this->ERROR_STATE = TOI_CORRUPTED;
		exit(1);
	}

	this->buffer.resize(filesize);
	const U64 readUntil = this->header.n_entries * sizeof(entry_type);

	if(readUntil % sizeof(entry_type) != 0){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "Mangled data!" << std::endl;
		this->ERROR_STATE = TOI_CORRUPTED;
		exit(1);
	}

	this->stream.read(this->buffer.data, readUntil);
	this->entries = reinterpret_cast<const entry_type*>(this->buffer.data);
	this->n_entries = readUntil / sizeof(entry_type);

	if(!(this->header.controller.sorted && this->header.controller.expanded)){
		if(this->stream.tellg() != filesize){
			std::cerr << Helpers::timestamp("ERROR", "TOI") << "Mangled data!" << std::endl;
			this->ERROR_STATE = TOI_CORRUPTED;
			exit(1);
		}
	} else {
		std::cerr << "is sorted and expanded" << std::endl;
		this->index = new index_type(this->header.n_contig, contigs);
		stream >> *this->index;

		if(!stream.good()){
			std::cerr << Helpers::timestamp("ERROR", "TOI") << "Corrupted data!" << std::endl;
			exit(1);
		}

		if(!stream.tellg() == filesize){
			std::cerr << Helpers::timestamp("ERROR", "TOI") << "Mangled data!" << std::endl;
			exit(1);
		}
	}

	this->ERROR_STATE = TOI_OK;
	return true;
}

bool TotempoleOutputReader::findOverlap(const U32 contigID, totempole_entry& intervals){
	if(this->ERROR_STATE != TOI_OK){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "No primary index available..." << std::endl;
		return false;
	}

	if(!this->header.controller.sorted){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "Index is not sorted..." << std::endl;
		return false;
	}

	if(this->index->getState() != TotempoleOutputSortedIndex::TOI_SORTED_ERROR::TOI_SORTED_OK){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "No sorted index available..." << std::endl;
		return false;
	}

	return(this->index->findOverlap(contigID, intervals));
}

bool TotempoleOutputReader::findOverlap(const U32 contigID, const U32 position, totempole_entry& intervals){
	if(this->ERROR_STATE != TOI_OK){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "No primary index available..." << std::endl;
		return false;
	}

	if(!this->header.controller.sorted){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "Index is not sorted..." << std::endl;
		return false;
	}

	if(this->index->getState() != TotempoleOutputSortedIndex::TOI_SORTED_ERROR::TOI_SORTED_OK){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "No sorted index available..." << std::endl;
		return false;
	}

	return(this->index->findOverlap(contigID, position, intervals));
}

bool TotempoleOutputReader::findOverlap(const U32 contigID, const U32 from, const U32 to, std::vector<totempole_entry>& intervals){
	if(this->ERROR_STATE != TOI_OK){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "No primary index available..." << std::endl;
		return false;
	}

	if(!this->header.controller.sorted){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "Index is not sorted..." << std::endl;
		return false;
	}

	if(this->index->getState() != TotempoleOutputSortedIndex::TOI_SORTED_ERROR::TOI_SORTED_OK){
		std::cerr << Helpers::timestamp("ERROR", "TOI") << "No sorted index available..." << std::endl;
		return false;
	}

	return(this->index->findOverlap(contigID, from, to, intervals));
}


} /* namespace Totempole */
} /* namespace Tomahawk */
