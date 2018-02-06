#ifndef TOMAHAWK_BASE_GENOTYPE_CONTAINER_RUNLENGTH_H_
#define TOMAHAWK_BASE_GENOTYPE_CONTAINER_RUNLENGTH_H_

#include <cassert>

#include "genotype_container_runlength_objects.h"
#include "genotype_objects.h"
#include "meta_entry.h"

namespace Tomahawk{
namespace Base{

template <class T>
class GenotypeContainerRunlength{
private:
    typedef GenotypeContainerRunlength           self_type;
    typedef GenotypeContainerRunlengthObjects<T> value_type;
    typedef value_type&                          reference;
    typedef const value_type&                    const_reference;
    typedef value_type*                          pointer;
    typedef const value_type*                    const_pointer;
    typedef std::ptrdiff_t                       difference_type;
    typedef std::size_t                          size_type;
	typedef MetaEntry<T>                 meta_type;

public:
	GenotypeContainerRunlength() :
		n_entries(0),
		__entries(nullptr)
	{}

	GenotypeContainerRunlength(const char* const genotype_buffer, const size_t l_buffer_length, const size_t n_entries, const meta_type* const meta_entries) :
		n_entries(n_entries),
		__entries(static_cast<pointer>(::operator new[](n_entries*sizeof(value_type))))
	{
		assert(n_entries > 0);
		assert(l_buffer_length % sizeof(T) == 0);

		size_t cumulative_position = 0;
		for(size_t i = 0; i < this->size(); ++i){
			new( &this->__entries[i] ) value_type( &genotype_buffer[cumulative_position], meta_entries[i].runs * sizeof(T) );
			cumulative_position += meta_entries[i].runs * sizeof(T);
			assert(this->__entries[i].size() > 0);
		}
		assert(cumulative_position == l_buffer_length);
	}

	~GenotypeContainerRunlength(){
		for(std::size_t i = 0; i < this->n_entries; ++i)
			((this->__entries + i)->~GenotypeContainerRunlengthObjects<T>)();

		::operator delete[](static_cast<void*>(this->__entries));
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
	inline reference at(const size_type& position){ return(this->__entries[position]); }
	inline const_reference at(const size_type& position) const{ return(this->__entries[position]); }
	inline reference operator[](const size_type& position){ return(this->__entries[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->__entries[position]); }
	inline pointer data(void){ return(this->__entries); }
	inline const_pointer data(void) const{ return(this->__entries); }
	inline reference front(void){ return(this->__entries[0]); }
	inline const_reference front(void) const{ return(this->__entries[0]); }
	inline reference back(void){ return(this->__entries[this->n_entries - 1]); }
	inline const_reference back(void) const{ return(this->__entries[this->n_entries - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }

	// Iterator
	inline iterator begin(){ return iterator(&this->__entries[0]); }
	inline iterator end()  { return iterator(&this->__entries[this->n_entries - 1]); }
	inline const_iterator begin()  const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->__entries[this->n_entries - 1]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->__entries[this->n_entries - 1]); }

private:
	size_type n_entries;
	pointer   __entries;
};

}
}



#endif /* TOMAHAWK_BASE_GENOTYPE_CONTAINER_RUNLENGTH_H_ */
