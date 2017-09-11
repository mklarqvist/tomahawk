
#ifndef SRC_TOTEMPOLE_TOTEMPOLEOUTPUTREADER_H_
#define SRC_TOTEMPOLE_TOTEMPOLEOUTPUTREADER_H_

#include "../support/MagicConstants.h"
#include "TotempoleMagic.h"
#include "../io/BasicBuffer.h"

namespace Tomahawk {
namespace Totempole {

class TotempoleOutputReader {
	typedef TotempoleOutputReader self_type;
	typedef Tomahawk::IO::TomahawkOutputSortHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> header_type;
	typedef TotempoleOutputEntry entry_type;
	typedef IO::BasicBuffer buffer_type;

	enum TOI_ERROR {TOI_OK, TOI_NO_EXIST, TOI_CORRUPTED};

public:
	TotempoleOutputReader() : n_entries(0), ERROR_STATE(TOI_OK), entries(nullptr){}
	~TotempoleOutputReader(){}

	bool Open(const std::string& input){
		this->stream = std::ifstream(input, std::ios::in | std::ios::binary | std::ios::ate);
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

	const entry_type& operator[](const U32 p) const{ return(this->entries[p]); }

	// Find data blocks mapping to these regions
	//bool findOverlap();

	bool getIsSorted(void) const{ return(this->header.controller.sorted); }
	const U32& size(void) const{ return(this->n_entries); }

public:
	TOI_ERROR ERROR_STATE;

private:
	U32 n_entries;
	std::ifstream stream;
	header_type header;
	buffer_type buffer;
	const entry_type* entries;
};

} /* namespace Totempole */
} /* namespace Tomahawk */

#endif /* SRC_TOTEMPOLE_TOTEMPOLEOUTPUTREADER_H_ */
