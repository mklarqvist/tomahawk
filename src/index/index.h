#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_

#include "index_entry.h"
#include "index_contig.h"
#include "index_container.h"
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
    	container_(data, l_data)
    {

    }

    inline const size_type& size(void) const{ return(this->container_.size()); }

    inline container_type& getContainer(void){ return(this->container_); }
	inline const container_type& getContainer(void) const{ return(this->container_); }

	inline void operator+=(const_reference entry){ this->container_ += entry; }

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& index){
		stream << index.getContainer();
		return(stream);
	}

private:
    container_type container_;
    // Todo: meta_index_container
};

}

#endif /* INDEX_INDEX_H_ */
