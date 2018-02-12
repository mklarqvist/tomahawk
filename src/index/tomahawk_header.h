#ifndef INDEX_TOMAHAWK_HEADER_H_
#define INDEX_TOMAHAWK_HEADER_H_

#include "index_header.h"
#include "index_entry.h"
#include "index_contig.h"
#include "index_container.h"
#include "../io/BasicBuffer.h"
#include "../algorithm/OpenHashTable.h"
#include "index.h"

namespace Tomahawk{

/**<
 * This container handles the header data for
 * a `twk` file
 */
class TomahawkHeader{
private:
	typedef TomahawkHeader             self_type;
	typedef Totempole::IndexHeader     header_type;
	typedef Totempole::IndexContigBase base_contig_type;
	typedef Totempole::IndexContig     contig_type;
    typedef Totempole::IndexEntry      entry_type;

    typedef std::ptrdiff_t             difference_type;
    typedef std::size_t                size_type;

    typedef IO::BasicBuffer            buffer_type;
    typedef Hash::OpenHashEntry<std::string, U32> hash_table;

public:

    // Accessors
    inline header_type& getHeader(void){ return(this->header_); }
    inline const header_type& getHeader(void) const{ return(this->header_); }
    inline std::string& getLiterals(void){ return(this->literals_); }
	inline const std::string& getLiterals(void) const{ return(this->literals_); }
	inline std::string& getSample(const U32 position){ return(this->samples_[position]); }
	inline const std::string& getSample(const U32 position) const{ return(this->samples_[position]); }

	inline const base_contig_type* getContigBase(const U32 contigID) const{ return(reinterpret_cast<const base_contig_type*>(&this->contigs_[contigID])); }

	// Updaters
	inline void addLiteral(const std::string& string){ this->literals_ += string; }

private:
	inline const bool validate(void) const{ return(this->header_.validate()); }
	bool BuildHashTables(void);

private:
    header_type    header_;
	std::string    literals_; // literal data
	contig_type*   contigs_;  // contig data
	std::string*   samples_;  // sample names
	hash_table*    contigsHashTable_; // contig name hash table
	hash_table*    sampleHashTable_;  // sample name hash table
};

}



#endif /* INDEX_TOMAHAWK_HEADER_H_ */
