#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_

#include "index_header.h"
#include "index_entry.h"
#include "index_contig.h"
#include "index_container.h"
#include "../algorithm/OpenHashTable.h"

namespace Tomahawk{

/**<
 * This container handles the index for `twk`
 */
class BlockIndex{
private:
	typedef BlockIndex                self_type;
	typedef Totempole::IndexHeader    header_type;
    typedef Totempole::IndexContig    contig_type;
    typedef Totempole::IndexEntry     entry_type;

    typedef std::ptrdiff_t            difference_type;
    typedef std::size_t               size_type;

    typedef IO::BasicBuffer           buffer_type;
    typedef Hash::OpenHashEntry<std::string, U32> hash_table;

public:



private:
    header_type    header;
	std::string    literals; // literal data
	contig_type*   contigs;  // contig data
	std::string*   samples;  // sample names
	entry_type*    entries;  // totempole entries data
	hash_table*    contigsHashTable; // contig name hash table
	hash_table*    sampleHashTable;  // sample name hash table
};

}



#endif /* INDEX_INDEX_H_ */
