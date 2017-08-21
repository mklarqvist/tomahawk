#ifndef BCFREADER_H_
#define BCFREADER_H_

#include <cassert>

#define BCF_ASSERT 1

#include "../BasicBuffer.h"
#include "../BGZFController.h"
#include "BCFEntry.h"

namespace Tomahawk {
namespace BCF {

class BCFReader{
	typedef BCFReader self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::BGZFController bgzf_controller_type;
	typedef IO::BGZFHeader bgzf_type;
	typedef VCF::VCFHeader header_type;
	typedef VCF::VCFHeaderContig contig_type;
	typedef BCFEntry entry_type;

public:
	BCFReader();
	~BCFReader();

	bool basicStats(const entry_type& entry);

	bool nextBlock(void);
	bool nextVariant(BCFEntry& entry);
	bool parseHeader(void);
	bool open(const std::string input);

public:
	std::ifstream stream;
	U64 filesize;
	U32 current_pointer;
	buffer_type buffer;
	buffer_type output_buffer;
	buffer_type header_buffer;
	bgzf_controller_type bgzf_controller;
	header_type header;
};

}
}

#endif /* BCFREADER_H_ */
