#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_

#include "index_contig.h"
#include "index_container.h"
#include "index_meta_container.h"
#include "../io/BasicBuffer.h"
#include "footer.h"

namespace Tomahawk{

/**<
 * This container handles the index entries for `twk` blocks: their
 * start and end IO positions and what genomic regions they cover.
 * The value type of this container are containers of entries.
 */
class Index{
private:
	typedef Index                     self_type;
	typedef Totempole::Footer         footer_type;
	typedef Totempole::IndexEntry     value_type;
	typedef Totempole::IndexContainer container_type;
	typedef Totempole::IndexMetaContainer meta_container_type;
	typedef Totempole::IndexMetaEntry meta_entry_type;
    typedef value_type&               reference;
    typedef const value_type&         const_reference;
    typedef value_type*               pointer;
    typedef const value_type*         const_pointer;
    typedef std::ptrdiff_t            difference_type;
    typedef std::size_t               size_type;
    typedef IO::BasicBuffer           buffer_type;

public:
    Index(){}
    ~Index(){}

    // Reading an index from a byte stream
    Index(const char* const data, const U32 l_data) :
    	meta_container_(data, *reinterpret_cast<const size_type* const>(data)*TWK_INDEX_META_ENTRY_SIZE+sizeof(size_type)),
    	container_(&data[this->meta_container_.size() * TWK_INDEX_META_ENTRY_SIZE + sizeof(size_type)], l_data - (this->meta_container_.size() * TWK_INDEX_META_ENTRY_SIZE + sizeof(size_type)))
    {

    }

    inline const size_type& size(void) const{ return(this->container_.size()); }

    inline container_type& getContainer(void){ return(this->container_); }
	inline const container_type& getContainer(void) const{ return(this->container_); }
	inline meta_container_type& getMetaContainer(void){ return(this->meta_container_); }
	inline const meta_container_type& getMetaContainer(void) const{ return(this->meta_container_); }


	inline void operator+=(const_reference entry){ this->container_ += entry; }

	bool buildMetaIndex(const U32 n_contigs){
		if(this->getContainer().size() == 0)
			return false;

		meta_entry_type reference_entry;
		reference_entry.index_begin = 0;
		reference_entry.index_end   = 1;
		reference_entry.min_position = this->getContainer()[0].min_position;
		reference_entry.max_position = this->getContainer()[0].max_position;
		reference_entry.n_variants = this->getContainer()[0].n_variants;
		reference_entry.uncompressed_size = this->getContainer()[0].uncompressed_size;
		U32 reference_contig = this->getContainer()[0].contigID;

		meta_entry_type* temp_entries = new meta_entry_type[n_contigs];

		for(U32 i = 1; i < this->getContainer().size(); ++i){
			if(this->getContainer()[i].contigID != reference_contig){
				if(this->getContainer()[i].contigID < reference_contig)
					continue;

				temp_entries[reference_contig] = reference_entry;
				reference_contig = this->getContainer()[i].contigID;
				reference_entry.index_begin = i;
				reference_entry.index_end = i + 1;
				reference_entry.min_position = this->getContainer()[i].min_position;
				reference_entry.max_position = this->getContainer()[i].max_position;
				reference_entry.n_variants = this->getContainer()[i].n_variants;
				reference_entry.uncompressed_size = this->getContainer()[i].uncompressed_size;

			} else {
				++reference_entry.index_end;
				reference_entry.max_position = this->getContainer()[i].max_position;
				reference_entry.n_variants += this->getContainer()[i].n_variants;
				reference_entry.uncompressed_size += this->getContainer()[i].uncompressed_size;
			}
		}
		temp_entries[reference_contig] = reference_entry;

		for(U32 i = 0; i < n_contigs; ++i){
			//std::cerr << temp_entries[i] << std::endl;
			this->meta_container_ += temp_entries[i];
		}
		delete [] temp_entries;

		return(true);
	}

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& index){
		stream << index.getMetaContainer();
		stream << index.getContainer();
		return(stream);
	}

private:
    meta_container_type meta_container_;
    container_type      container_;
};

}

#endif /* INDEX_INDEX_H_ */
