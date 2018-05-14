#ifndef BCFREADER_H_
#define BCFREADER_H_

#include <cassert>

#define BCF_ASSERT 1

#include "../basic_buffer.h"
#include "../compression/bgzf_controller.h"
#include "bcf_entry.h"
#include "../vcf/vcf_header.h"

namespace tomahawk {
namespace bcf {

class BCFReader{
	typedef BCFReader            self_type;
	typedef io::BasicBuffer      buffer_type;
	typedef io::BGZFController   bgzf_controller_type;
	typedef io::BGZFHeader       bgzf_type;
	typedef vcf::VCFHeader       header_type;
	typedef vcf::VCFHeaderContig contig_type;
	typedef BCFEntry             entry_type;

public:
	BCFReader();
	~BCFReader();

	bool basicStats(const entry_type& entry);

	bool nextBlock(void);
	bool nextVariant(BCFEntry& entry);
	bool parseHeader(void);
	bool open(const std::string input);

public:
	std::ifstream        stream;
	U64                  filesize;
	U32                  current_pointer;
	buffer_type          buffer;
	buffer_type          header_buffer;
	bgzf_controller_type bgzf_controller;
	header_type          header;
};

}
}

#endif /* BCFREADER_H_ */
