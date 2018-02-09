#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_

#include "index_header.h"
#include "index_entry.h"
#include "index_contig.h"
#include "index_container.h"

namespace Tomahawk{

class Index{
private:
	typedef Index self_type;
	typedef Totempole::IndexHeader    header_type;
	typedef Totempole::IndexContig    contig_type;
	typedef Totempole::IndexEntry     entry_type;
	typedef Totempole::IndexContainer value_type;
    typedef value_type&               reference;
    typedef const value_type&         const_reference;
    typedef value_type*               pointer;
    typedef const value_type*         const_pointer;
    typedef std::ptrdiff_t            difference_type;
    typedef std::size_t               size_type;
    typedef IO::BasicBuffer           buffer_type;

public:

private:
    size_type    n_entries;
    header_type  header;
    contig_type* contigs;
    pointer      __entries;
};

}



#endif /* INDEX_INDEX_H_ */
