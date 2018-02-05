#ifndef TOMAHAWK_BASE_GENOTYPE_CONTAINER_REFERENCE_H_
#define TOMAHAWK_BASE_GENOTYPE_CONTAINER_REFERENCE_H_

#include <cstring>  // size_t, ptrdiff_t

#include "../../support/TypeDefinitions.h"
#include "../../totempole/TotempoleEntry.h"
#include "genotype_container_bitvector.h"

namespace Tomahawk{
namespace Base{

/**<
 * Special genotype container for both run-length encoded
 * and bit-encoded genotypes using unaligned memory directly
 * interpreted from type-casts at compile time and has some
 * psuedo-iterator capabilities.
 * This works because there is no random access in the pairwise
 * comparator functions. Upper triagonal comparisons can be done
 * by strictly iterating forward.
 *
 * Advanced use only: requires down-casting into the template
 * primitive type that the given Tomahawk file has been encoded
 * in.
 * Data can only be interpreted invoking the standard constructor
 * All other constructors have non-standard meaning.
 * 1) Copy constructor: copies the pointer addresses and iterator
 *                      positions only.
 * 2) Assignment operator: copies the iterator position only
 *
 * For example:
 * GenotypeContainerReference gt;
 * GenotypeContainerReferenceImpl<U16>& impl = *reinterpret_cast<GenotypeContainerReferenceImpl<U16>*>(&gt);
 */
class GenotypeContainerReference{
private:
	typedef GenotypeContainerReference  self_type;

protected:
	typedef Totempole::TotempoleEntry   header_entry;
	typedef GenotypeContainerBitvector  container_bitvector_type;
	typedef Base::GenotypeBitvector<>   genotype_bitvector_type;
	typedef std::size_t                 size_type;

public:
	GenotypeContainerReference();
	virtual ~GenotypeContainerReference();

	//
	void operator=(const self_type& other); // copy pointers only!

	// Accessor
	inline const header_entry& getTotempole(void) const{ return(this->index_entry); }
	inline const genotype_bitvector_type& getBitvector(const U32& position) const{ return(this->bit_vectors->at(position)); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }
	inline void resetIterator(void){ this->iterator_position = 0; }

	// Accessor
	inline const char* const data(void) const{ return(this->genotype_data); }

	//
	virtual void operator++(void) =0;
	virtual void operator--(void) =0;

protected:
	size_type                 n_entries;
	size_type                 iterator_position;
	const char*               genotype_data;
	const char*               meta_data;
	header_entry*             index_entry;
	container_bitvector_type* bit_vectors;
};

/**<
 * Template T represents the primitive that genotypes are
 * encoded in.
 */
template <class T>
class GenotypeContainerReferenceImpl : public GenotypeContainerReference{
private:
	typedef GenotypeContainerReferenceImpl       self_type;
	typedef GenotypeContainerRunlengthObjects<T> value_type;
	typedef value_type&                          reference;
	typedef const value_type&                    const_reference;
	typedef value_type*                          pointer;
	typedef const value_type*                    const_pointer;
	typedef std::ptrdiff_t                       difference_type;
	typedef std::size_t                          size_type;
	typedef TomahawkEntryMeta<T>                 meta_type;

public:
	GenotypeContainerReferenceImpl();
	~GenotypeContainerReferenceImpl();

	inline const meta_type& getMeta(const U32& position) const{ return(*reinterpret_cast<const meta_type* const>(&this->meta_data[position * (TOMAHAWK_ENTRY_META_SIZE + sizeof(T))])); }
	inline const meta_type& currentMeta(void) const{ return(*reinterpret_cast<const meta_type* const>(&this->meta_data[this->iterator_position * (TOMAHAWK_ENTRY_META_SIZE + sizeof(T))])); }
	inline const_pointer current(void) const{ return(reinterpret_cast<const_pointer>(&this->genotype_data[this->iterator_position * sizeof(value_type)])); }

	// Psuedo-iterator functionality
	inline void operator++(void){ this->iterator_position += this->currentMeta().runs; }
	inline void operator--(void){ this->iterator_position -= this->currentMeta().runs; }
};

}
}



#endif /* TOMAHAWK_BASE_GENOTYPE_CONTAINER_REFERENCE_H_ */
