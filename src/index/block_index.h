#ifndef INDEX_BLOCK_INDEX_H_
#define INDEX_BLOCK_INDEX_H_

#include "index_header.h"
#include "index_entry.h"
#include "index_contig.h"
#include "index_container.h"

namespace Tomahawk{

/**<
 * This container handles the index entries for `twk` blocks: their
 * start and end IO positions and what genomic regions they cover.
 * The value type of this container are containers of entries.
 */
class BlockIndex{
private:
	typedef BlockIndex                self_type;
	typedef Totempole::IndexHeader    header_type;
	typedef Totempole::IndexEntry     value_type;
	typedef Totempole::IndexContainer container_type;
    typedef value_type&               reference;
    typedef const value_type&         const_reference;
    typedef value_type*               pointer;
    typedef const value_type*         const_pointer;
    typedef std::ptrdiff_t            difference_type;
    typedef std::size_t               size_type;
    typedef IO::BasicBuffer           buffer_type;

public:

private:
    header_type    header;
    container_type container;
};

}

#endif /* INDEX_BLOCK_INDEX_H_ */
