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
	typedef Index                         self_type;
	typedef Totempole::Footer             footer_type;
	typedef Totempole::IndexEntry         value_type;
	typedef Totempole::IndexContainer     container_type;
	typedef Totempole::IndexMetaContainer meta_container_type;
	typedef Totempole::IndexMetaEntry     meta_entry_type;
    typedef value_type&                   reference;
    typedef const value_type&             const_reference;
    typedef value_type*                   pointer;
    typedef const value_type*             const_pointer;
    typedef std::ptrdiff_t                difference_type;
    typedef std::size_t                   size_type;
    typedef IO::BasicBuffer               buffer_type;

public:
    Index();
    ~Index();

    // Reading an index from a byte stream
    Index(const char* const data, const U32 l_data);

    // Capacity
    inline const size_type& size(void) const{ return(this->container_.size()); }
    inline const size_type& sizeMeta(void) const{ return(this->meta_container_.size()); }

    // Accessors
    inline container_type& getContainer(void){ return(this->container_); }
	inline const container_type& getContainer(void) const{ return(this->container_); }
	inline meta_container_type& getMetaContainer(void){ return(this->meta_container_); }
	inline const meta_container_type& getMetaContainer(void) const{ return(this->meta_container_); }

	// Overloaded
	inline void operator<<(const_reference entry){ this->container_ += entry; }
	inline void operator+=(const_reference entry){ this->container_ += entry; }

	/**<
	 *
	 * @param n_contigs
	 * @return
	 */
	bool buildMetaIndex(const U32 n_contigs);

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
