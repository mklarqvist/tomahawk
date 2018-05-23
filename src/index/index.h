#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_

#include "io/basic_buffer.h"
#include "index_contig.h"
#include "index_container.h"
#include "index_meta_container.h"
#include "footer.h"

namespace tomahawk{

/**<
 * Index controller for bit flags
 */
struct IndexController{
public:
	typedef IndexController self_type;

public:
	IndexController() :
		isSorted(false),
		isPartialSorted(false),
		unused(0)
	{}

	IndexController(const char* const data){ memcpy(this, data, sizeof(BYTE)); }

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& controller){
		stream.write((const char*)&controller, sizeof(BYTE));
		return stream;
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& controller){
		stream.read((char*)&controller, sizeof(BYTE));
		return stream;
	}

public:
	BYTE isSorted:        1,
         isPartialSorted: 1,
	     unused:          6;
};

/**<
 * This container handles the index entries for `twk` blocks: their
 * start and end IO positions and what genomic regions they cover.
 * The value type of this container are containers of entries.
 */
class Index{
private:
	typedef Index                         self_type;
	typedef totempole::Footer             footer_type;
	typedef totempole::IndexEntry         value_type;
	typedef totempole::IndexContainer     container_type;
	typedef totempole::IndexMetaContainer meta_container_type;
	typedef totempole::IndexMetaEntry     meta_entry_type;
    typedef value_type&                   reference;
    typedef const value_type&             const_reference;
    typedef value_type*                   pointer;
    typedef const value_type*             const_pointer;
    typedef std::ptrdiff_t                difference_type;
    typedef std::size_t                   size_type;
    typedef io::BasicBuffer               buffer_type;
    typedef IndexController               controller_type;

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
	inline controller_type& getController(void){ return(this->controller_); }
	inline const controller_type& getController(void) const{ return(this->controller_); }

	// Setters
	inline void setSorted(const bool yes){ this->controller_.isSorted = yes; }
	inline void setPartialSorted(const bool yes){ this->controller_.isPartialSorted = yes; }

	// Getters
	inline const bool isSorted(void) const{ return(this->controller_.isSorted); }
	inline const bool isPartialSorted(void) const{ return(this->controller_.isPartialSorted); }
	inline const U64  totalBytes(void) const{
		U64 total_bytes = 0;
		for(size_t i = 0; i < this->getContainer().size(); ++i)
			total_bytes += this->getContainer().at(i).sizeBytes();

		return(total_bytes);
	}

	// Overloaded
	inline void operator<<(const_reference entry){ this->container_ += entry; }
	inline void operator+=(const_reference entry){ this->container_ += entry; }

	/**<
	 * Constructs the index of index if the data is sorted
	 * @param n_contigs Number of contigs in the file
	 * @return          Returns TRUE upon success or FALSE otherwise
	 */
	bool buildMetaIndex(const U32 n_contigs);

	/**<
	 *
	 * @param contigID
	 * @param fromPos
	 * @param toPos
	 * @return
	 */
	std::vector<U32> findOverlaps(const U32 contigID, const U32 fromPos, const U32 toPos) const{
		if(contigID > this->sizeMeta())
			return std::vector<U32>();

		const U32 blockFrom = this->meta_container_[contigID].index_begin;
		const U32 blockTo   = this->meta_container_[contigID].index_end;
		std::vector<U32> ret;

		for(U32 i = blockFrom; i < blockTo; ++i){
			// [a, b] overlaps with [x, y] iff b > x and a < y.
			// a = fromPos
			// b = toPos
			// x = this->container_[i].min_position
			// y = this->container_[i].max_position
			if(toPos >= this->container_[i].min_position && fromPos < this->container_[i].max_position)
				ret.push_back(i);

			// Cannot extend any more
			if(this->container_[i].min_position > toPos) break;
		}

		return ret;
	}

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& index){
		stream << index.getController();
		stream << index.getMetaContainer();
		stream << index.getContainer();
		return(stream);
	}

private:
	controller_type     controller_;
    meta_container_type meta_container_;
    container_type      container_;
};

}

#endif /* INDEX_INDEX_H_ */
