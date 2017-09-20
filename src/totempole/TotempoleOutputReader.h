
#ifndef SRC_TOTEMPOLE_TOTEMPOLEOUTPUTREADER_H_
#define SRC_TOTEMPOLE_TOTEMPOLEOUTPUTREADER_H_

#include "../support/MagicConstants.h"
#include "TotempoleMagic.h"
#include "../io/BasicBuffer.h"
#include "TotempoleOutputSortedIndex.h"

namespace Tomahawk {
namespace Totempole {

class TotempoleOutputReader {
	typedef TotempoleOutputReader self_type;
	typedef Tomahawk::IO::TomahawkOutputSortHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> header_type;
	typedef TotempoleOutputEntry entry_type;
	typedef IO::BasicBuffer buffer_type;

	enum TOI_ERROR {TOI_OK, TOI_NO_EXIST, TOI_CORRUPTED};

public:
	TotempoleOutputReader();
	~TotempoleOutputReader();

	bool Open(const std::string& input);

	inline const entry_type& operator[](const U32 p) const{ return(this->entries[p]); }
	inline bool getIsSorted(void) const{ return(this->header.controller.sorted); }
	inline const U32& size(void) const{ return(this->n_entries); }

	// Find data blocks mapping to these regions
	bool findOverlap(const S32 contigID);
	bool findOverlap(const S32 contigID, const U32 position);
	bool findOverlap(const S32 contigID, const U32 from, const U32 to);


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
