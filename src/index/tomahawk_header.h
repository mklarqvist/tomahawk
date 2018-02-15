#ifndef INDEX_TOMAHAWK_HEADER_H_
#define INDEX_TOMAHAWK_HEADER_H_

#include <cassert>

#include "../algorithm/open_hashtable.h"
#include "index_contig.h"
#include "../io/BasicBuffer.h"
#include "../tomahawk/tomahawk_magic_header.h"
#include "../io/compression/TGZFController.h"

namespace Tomahawk{

/**<
 * This container handles the header data for
 * a `twk` file
 */
class TomahawkHeader{
public:
	typedef TomahawkHeader                    self_type;
	typedef Totempole::HeaderContig           contig_type;
    typedef Base::TomahawkMagicHeader         magic_type;
    typedef IO::BasicBuffer                   buffer_type;
    typedef Hash::HashTable<std::string, S32> hash_table;
    typedef IO::TGZFController                compressor_type;

public:
    TomahawkHeader(void);

    // Standard dtor
    ~TomahawkHeader(void);

    // Open and close functions
    int open(std::istream& stream = std::cin);
    int write(std::ostream& stream = std::cout);

    // Accessors
    inline std::string& getLiterals(void){ return(this->literals_); }
	inline const std::string& getLiterals(void) const{ return(this->literals_); }
	inline std::string& getSample(const U32 position){ return(this->sample_names_[position]); }
	inline const std::string& getSample(const U32 position) const{ return(this->sample_names_[position]); }

	const bool getSample(const std::string& sample_name, const std::string*& return_target) const;
	const bool getContigName(const std::string& contig_name, const std::string*& return_target) const;
	const bool getContig(const std::string& contig_name, const contig_type*& return_target) const;
	const S32 getContigID(const std::string& contig_name) const;

	inline magic_type& getMagic(void){ return(this->magic_); }
	inline const magic_type& getMagic(void) const{ return(this->magic_); }

	// Updater
	inline void addLiteral(const std::string& string){ this->literals_ += string; }
	inline const bool validate(void) const{ return(this->magic_.validate()); }

private:
	bool BuildHashTables(void);
	const U32 DetermineUncompressedSize(void) const;

public:
	magic_type     magic_;              // magic header
	std::string    literals_;           // literal data
	contig_type*   contigs_;            // contig data
	std::string*   sample_names_;       // sample names
	hash_table*    contigs_hash_table_; // contig name hash table
	hash_table*    sample_hash_table_;  // sample name hash table
};

}

#endif /* INDEX_TOMAHAWK_HEADER_H_ */
