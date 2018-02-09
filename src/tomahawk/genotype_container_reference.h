#ifndef TOMAHAWK_BASE_GENOTYPE_CONTAINER_REFERENCE_H_
#define TOMAHAWK_BASE_GENOTYPE_CONTAINER_REFERENCE_H_

#include <cstring>  // size_t, ptrdiff_t

#include "../support/type_definitions.h"
#include "../index/index_entry.h"
#include "genotype_container_bitvector.h"

namespace Tomahawk{
namespace Base{

/**<
 * Special genotype container for both run-length encoded
 * and bit-encoded genotypes using unaligned memory directly
 * interpreted from type-casts at compile time and has some
 * psuedo-iterator capabilities.
 * This works because there is no random access in the pairwise
 * comparator functions. Upper triangular comparisons can be done
 * by strictly iterating forward.
 *
 * Data can only be interpreted by invoking the standard constructor.
 * All other constructors have non-standard meaning.
 * 1) Copy constructor: copies the pointer addresses and iterator
 *                      positions only.
 * 2) Assignment operator: copies the iterator position only
 */
template <class T>
class GenotypeContainerReference{
private:
	typedef GenotypeContainerReference  self_type;

protected:
	typedef Totempole::IndexEntry          header_entry_type;
	typedef GenotypeContainerBitvector     container_bitvector_type;
	typedef Base::GenotypeBitvector<>      genotype_bitvector_type;
	typedef Support::GenotypeDiploidRun<T> value_type;
	typedef value_type&                    reference;
	typedef const value_type&              const_reference;
	typedef value_type*                    pointer;
	typedef const value_type*              const_pointer;
	typedef std::ptrdiff_t                 difference_type;
	typedef std::size_t                    size_type;
	typedef MetaEntry<T>                   meta_type;

public:
	GenotypeContainerReference() :
		n_entries(0),
		iterator_position_meta(0),
		iterator_position_runs(0),
		owns_bitvectors(true),
		meta_entries(nullptr),
		genotype_entries(nullptr),
		index_entry(nullptr),
		bit_vectors(nullptr)
	{
	}

	GenotypeContainerReference(const char* const data, const size_t l_data, const header_entry_type& index_entry, const size_t n_samples) :
		n_entries(index_entry.n_variants),
		iterator_position_meta(0),
		iterator_position_runs(0),
		owns_bitvectors(true),
		meta_entries(reinterpret_cast<const meta_type* const>(data)),
		genotype_entries(reinterpret_cast<const_pointer>(&data[this->n_entries * (TOMAHAWK_ENTRY_META_SIZE + sizeof(T))])),
		index_entry(&index_entry),
		bit_vectors(nullptr)
	{
		this->bit_vectors = new container_bitvector_type();
		this->bit_vectors->Build(this->genotype_entries, this->meta_entries, this->size(), n_samples);
	}

	// Copy ctor: copies iterator positions and pointers
	GenotypeContainerReference(const self_type& other) :
		n_entries(other.n_entries),
		iterator_position_meta(other.iterator_position_meta),
		iterator_position_runs(other.iterator_position_runs),
		owns_bitvectors(false),
		meta_entries(other.meta_entries),
		genotype_entries(other.genotype_entries),
		index_entry(other.index_entry),
		bit_vectors(other.bit_vectors)
	{

	}

	~GenotypeContainerReference(){
		if(this->owns_bitvectors)
			delete this->bit_vectors;
	}

	// // copy pointers only!
	void operator=(const self_type& other){
		this->iterator_position_runs = other.iterator_position_runs;
		this->iterator_position_meta = other.iterator_position_meta;
	}

	// Accessor
	inline const header_entry_type& getTotempole(void) const{ return(*this->index_entry); }
	inline const genotype_bitvector_type& getBitvector(const U32& position) const{ return(this->bit_vectors->at(position)); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }
	inline void resetIterator(void){ this->iterator_position_runs = 0; this->iterator_position_meta = 0; }

	// Accessor
	inline const char* const data(void) const{ return(this->genotype_entries); }
	inline const char* const meta_data(void) const{ return(this->meta_entries); }

	inline const meta_type& getMeta(const U32& position) const{ return(this->meta_entries[position]); }
	inline const meta_type& currentMeta(void) const{ return(this->meta_entries[this->iterator_position_meta]); }
	inline const_pointer current(void) const{ return(&this->genotype_entries[this->iterator_position_runs]); }
	inline const genotype_bitvector_type& currentBitvector(void) const{ return(this->bit_vectors->at(this->iterator_position_meta)); }


	// Psuedo-iterator functionality
	inline void operator++(void){ this->iterator_position_runs += this->currentMeta().runs; ++this->iterator_position_meta; }
	inline void operator--(void){ this->iterator_position_runs -= this->currentMeta().runs; ++this->iterator_position_meta; }
	inline const_reference operator[](const U32& position) const{ return(this->genotype_entries[this->iterator_position_runs + position]); }
	inline const_reference at(const U32& position) const{ return(this->genotype_entries[this->iterator_position_runs + position]); }

protected:
	size_type                 n_entries;
	size_type                 iterator_position_meta;
	size_type                 iterator_position_runs;
	bool                      owns_bitvectors;
	const meta_type*          meta_entries;
	const_pointer             genotype_entries;
	const header_entry_type*  index_entry;
	container_bitvector_type* bit_vectors;
};

}
}

#endif /* TOMAHAWK_BASE_GENOTYPE_CONTAINER_REFERENCE_H_ */
