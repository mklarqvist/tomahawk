#ifndef TOMAHAWK_BASE_GENOTYPE_CONTAINER_H_
#define TOMAHAWK_BASE_GENOTYPE_CONTAINER_H_

#include <cstring>  // size_t, ptrdiff_t
#include <iterator> // forward_iterator_tag

#include "support/type_definitions.h"
#include "index/index_entry.h"
#include "genotype_container_bitvector.h"
#include "genotype_container_runlength.h"

namespace tomahawk{
namespace base{

template <class T>
class GenotypeContainer{
private:
	typedef GenotypeContainer                    self_type;
	typedef GenotypeContainerBitvector           container_bitvector_type;
	typedef GenotypeContainerRunlength<T>        container_runlength_type;
	typedef GenotypeContainerRunlengthObjects<T> genotype_runlength_type;
	typedef base::GenotypeBitvector<>            genotype_bitvector_type;
	typedef MetaEntry                            meta_type;
	typedef totempole::IndexEntry                header_entry;
    typedef genotype_runlength_type              value_type;
    typedef value_type&                          reference;
    typedef const value_type&                    const_reference;
    typedef value_type*                          pointer;
    typedef const value_type*                    const_pointer;
    typedef std::ptrdiff_t                       difference_type;
    typedef std::size_t                          size_type;

public:
	GenotypeContainer(const char* const data_buffer, const size_t l_buffer_length, const header_entry& support, const U64& n_samples) :
		n_entries(support.n_variants),
		index_entry(support), // invoke copy ctor
		meta_entries(static_cast<meta_type*>(::operator new[](this->n_entries*sizeof(meta_type)))),
		container_runlength(nullptr),
		container_bitvector(nullptr)
	{
		if(l_buffer_length == 0)
			return;

		// Interpret meta entries
		size_t cumulative_position = 0;
		size_t genotype_cost = 0;

		for(size_t i = 0; i < this->size(); ++i){
			new( &this->meta_entries[i] ) meta_type( &data_buffer[cumulative_position] );
			cumulative_position += TOMAHAWK_ENTRY_META_SIZE;
			genotype_cost += meta_entries[i].runs*sizeof(T);
		}
		assert(cumulative_position + genotype_cost == l_buffer_length);

		// Interpret run lengths
		this->container_runlength = new container_runlength_type(&data_buffer[cumulative_position], l_buffer_length - cumulative_position, this->size(), this->meta_entries);

		// Interpret bit vectors
		this->container_bitvector = new container_bitvector_type(*this->container_runlength, n_samples);
	}

	~GenotypeContainer(){
		for(size_type i = 0; i < this->size(); ++i)
			((this->meta_entries + i)->~MetaEntry)();

		::operator delete[](static_cast<void*>(this->meta_entries));
		delete this->container_runlength;
		delete this->container_bitvector;
	}

	class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const{ return *ptr_; }
		pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	class const_iterator{
	private:
		typedef const_iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		const_reference operator*() const{ return *ptr_; }
		const_pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	// Element access
	inline reference at(const size_type& position){ return(this->container_runlength->at(position)); }
	inline const_reference at(const size_type& position) const{ return(this->container_runlength->at(position)); }
	inline reference operator[](const size_type& position){ return(this->container_runlength->at(position)); }
	inline const_reference operator[](const size_type& position) const{ return(this->container_runlength->at(position)); }
	inline pointer data(void){ return(this->container_runlength); }
	inline const_pointer data(void) const{ return(this->container_runlength); }
	inline reference front(void){ return(this->container_runlength->at(0)); }
	inline const_reference front(void) const{ return(this->container_runlength->at(0)); }
	inline reference back(void){ return(this->container_runlength->at(this->n_entries - 1)); }
	inline const_reference back(void) const{ return(this->container_runlength->at(this->n_entries - 1)); }

	// Accessor
	inline const meta_type& getMeta(const U32& position) const{ return(this->meta_entries[position]); }
	inline const header_entry& getTotempole(void) const{ return(this->index_entry); }
	inline const genotype_bitvector_type& getBitvector(const U32& position) const{ return(this->container_bitvector->at(position)); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }

	// Iterator
	inline iterator begin(){ return iterator(&this->container_runlength[0]); }
	inline iterator end()  { return iterator(&this->container_runlength[this->n_entries - 1]); }
	inline const_iterator begin()  const{ return const_iterator(&this->container_runlength[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->container_runlength[this->n_entries - 1]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->container_runlength[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->container_runlength[this->n_entries - 1]); }

private:
	size_type                 n_entries;
	header_entry              index_entry;
	meta_type*                meta_entries;
	container_runlength_type* container_runlength;
	container_bitvector_type* container_bitvector;
};

}
}


#endif /* TOMAHAWK_BASE_GENOTYPE_CONTAINER_H_ */
