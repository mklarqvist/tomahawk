#ifndef TOMAHAWK_BASE_TWK_READER_IMPLEMENTATION_H_
#define TOMAHAWK_BASE_TWK_READER_IMPLEMENTATION_H_

#include "genotype_container.h"

namespace tomahawk{

template <class T>
class TomahawkReaderImpl{
private:
	typedef TomahawkReaderImpl         self_type;
	typedef base::GenotypeContainer<T> value_type;
    typedef value_type&                reference;
    typedef const value_type&          const_reference;
    typedef value_type*                pointer;
    typedef const value_type*          const_pointer;
    typedef std::ptrdiff_t             difference_type;
    typedef std::size_t                size_type;
	typedef MetaEntry                  meta_type;
	typedef totempole::IndexEntry      header_entry;
	typedef totempole::IndexEntry      support_type;

public:
	TomahawkReaderImpl(const U64 n_samples) :
		n_entries(0),
		n_capacity(0),
		n_samples(n_samples),
		__entries(nullptr)
	{

	}

	TomahawkReaderImpl(const U64 n_samples, const size_t n_capacity) :
		n_entries(0),
		n_capacity(n_capacity),
		n_samples(n_samples),
		__entries(static_cast<pointer>(::operator new[](this->n_capacity*sizeof(value_type))))
	{

	}

	~TomahawkReaderImpl(){
		for(size_type i = 0; i < this->size(); ++i)
			((this->__entries + i)->~value_type)();

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
	inline const size_type& capacity(void) const{ return(this->n_capacity); }

	// Iterator
	inline iterator begin(){ return iterator(&this->__entries[0]); }
	inline iterator end()  { return iterator(&this->__entries[this->n_entries - 1]); }
	inline const_iterator begin()  const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->__entries[this->n_entries - 1]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->__entries[this->n_entries - 1]); }

	/**<
	 * Add a `TWK` block as a reference container to the meta container
	 * @param data    Input data reference
	 * @param l_data  Length of input data
	 * @param support Paired `TWI` header with this block
	 * @return        Returns TRUE upon success or FALSE otherwise
	 */
	bool addDataBlock(const char* const data, const size_t l_data, const support_type& support){
		// Container is full
		// Resize is required
		if(this->n_entries + 1 == this->n_capacity || this->capacity() == 0)
			return false;

		//std::cerr << "constructing new @ " << this->n_entries << "/" << this->n_capacity << " and samples: " << this->n_samples << std::endl;
		new( &this->__entries[this->n_entries] ) value_type( data, l_data, support, this->n_samples );
		++this->n_entries;
		return true;
	}

	/**<
	 * Counts the number of variants referenced in this meta container
	 * @return Returns the total number of variants in all containers
	 */
	const U64 countVariants(void) const{
		U64 n_total = 0;
		for(U32 i = 0; i < this->size(); ++i)
			n_total += this->at(i).getTotempole().size();

		return(n_total);
	}

	const U64& numberSamples(void) const{ return(this->n_samples); }

private:
	size_type n_entries;
	size_type n_capacity;
	U64       n_samples;
	pointer   __entries;
};

}

#endif /* TOMAHAWK_BASE_TWK_READER_IMPLEMENTATION_H_ */
