#ifndef TOMAHAWK_BASE_GENOTYPE_CONTAINER_RUNLENGTH_OBJECTS_H_
#define TOMAHAWK_BASE_GENOTYPE_CONTAINER_RUNLENGTH_OBJECTS_H_

#include "genotype_objects.h"

namespace Tomahawk{
namespace Base{

/**<
 * Primary run-length encoded genotype objects used in Tomahawk.
 * These higher-order objects are encoded in fixed-width basic
 * unsigned primitives.
 * The template parameter represents the encoded primitive type
 * for a given Tomahawk file.
 */
template <class T>
class GenotypeContainerRunlengthObjects{
private:
    typedef GenotypeContainerRunlengthObjects self_type;
    typedef Support::GenotypeDiploidRun<T>    value_type;
    typedef value_type&                       reference;
    typedef const value_type&                 const_reference;
    typedef value_type*                       pointer;
    typedef const value_type*                 const_pointer;
    typedef std::ptrdiff_t                    difference_type;
    typedef std::size_t                       size_type;

public:
    GenotypeContainerRunlengthObjects() :
    	n_entries(0),
		__entries(nullptr)
	{}

    GenotypeContainerRunlengthObjects(const char* const genotype_buffer, const size_t l_buffer_length) :
    	n_entries(l_buffer_length / sizeof(T)),
		__entries(new value_type[this->n_entries])
    {
    	assert(l_buffer_length % sizeof(T) == 0);
    	memcpy(this->__entries, genotype_buffer, l_buffer_length);
    }

    ~GenotypeContainerRunlengthObjects(){
    	delete [] this->__entries;
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



#endif /* TOMAHAWK_BASE_GENOTYPE_CONTAINER_RUNLENGTH_OBJECTS_H_ */
