
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
	typedef Tomahawk::IO::TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> header_type;
	typedef TotempoleOutputEntry entry_type;
	typedef IO::BasicBuffer buffer_type;
	typedef TotempoleOutputSortedIndex index_type;
	typedef Totempole::TotempoleContigBase contig_type;
	typedef TotempoleOutputSortedEntry totempole_entry;

public:
	enum TOI_ERROR {TOI_OK, TOI_NO_EXIST, TOI_CORRUPTED, TOI_INIT};

public:
	TotempoleOutputReader();
	~TotempoleOutputReader();

	bool Open(const std::string& input, const contig_type* contigs);

	inline const entry_type& operator[](const U32 p) const{ return(this->entries[p]); }
	inline bool getIsSorted(void) const{ return(this->header.controller.sorted); }
	inline bool getIsSortedExpanded(void) const{ return(this->header.controller.sorted && this->header.controller.expanded); }
	inline const U32& size(void) const{ return(this->n_entries); }

	// TOI dispatchers
	// Find data blocks mapping to these regions
	bool findOverlap(const U32 contigID, totempole_entry& intervals);
	bool findOverlap(const U32 contigID, const U32 position, totempole_entry& intervals);
	bool findOverlap(const U32 contigID, const U32 from, const U32 to, std::vector<totempole_entry>& intervals);

	header_type& getHeader(void){ return(this->header); }
	const index_type* const getIndex(void) const{ return(this->index); }

private:
	bool __stateCheck(void) const;

public:
	TOI_ERROR ERROR_STATE;

private:
	U32 n_entries;
	std::ifstream stream;
	header_type header;
	buffer_type buffer;
	const entry_type* entries;
	index_type* index;
};

} /* namespace Totempole */
} /* namespace Tomahawk */

#endif /* SRC_TOTEMPOLE_TOTEMPOLEOUTPUTREADER_H_ */
