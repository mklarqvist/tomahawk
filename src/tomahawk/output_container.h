#ifndef TOMAHAWK_BASE_OUTPUT_CONTAINER_H_
#define TOMAHAWK_BASE_OUTPUT_CONTAINER_H_

#include <cassert>

#include "two/output_entry.h"

namespace Tomahawk{

class OutputContainer{
private:
    typedef IO::OutputEntry    value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;
    typedef std::ptrdiff_t     difference_type;
    typedef std::size_t        size_type;
    typedef IO::BasicBuffer    buffer_type;

public:
    OutputContainer() :
    	n_entries(0),
		n_capacity(0),
		__entries(nullptr)
	{

	}

    OutputContainer(const size_t capacity) :
    	n_entries(0),
		n_capacity(capacity),
		__entries(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
    {

    }

    OutputContainer(char* const data, const U64 l_data) :
    	n_entries(l_data / sizeof(value_type)),
		n_capacity(n_entries),
		__entries(static_cast<pointer>(::operator new[](this->size()*sizeof(value_type))))
	{
		assert(l_data % sizeof(value_type) == 0);

		U32 cumulative_position = 0;
		for(size_t i = 0; i < this->size(); ++i){
			new( &this->__entries[i] ) value_type( &data[cumulative_position] );
			cumulative_position += sizeof(value_type);
		}
		assert(cumulative_position == l_data);
	}

    OutputContainer(const buffer_type& data_buffer) :
		n_entries(data_buffer.size() / sizeof(value_type)),
		n_capacity(n_entries),
		__entries(static_cast<pointer>(::operator new[](this->size()*sizeof(value_type))))
	{
		assert(data_buffer.size() % sizeof(value_type) == 0);

		U32 cumulative_position = 0;
		for(size_t i = 0; i < this->size(); ++i){
			new( &this->__entries[i] ) value_type( &data_buffer[cumulative_position] );
			cumulative_position += sizeof(value_type);
		}
		assert(cumulative_position == data_buffer.size());
	}

    ~OutputContainer(){
    	for(size_type i = 0; i < this->size(); ++i)
			((this->__entries + i)->~OutputEntry)();

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
	inline const size_type& capacity(void) const{ return(this->capacity()); }

	// Iterator
	inline iterator begin(){ return iterator(&this->__entries[0]); }
	inline iterator end()  { return iterator(&this->__entries[this->n_entries - 1]); }
	inline const_iterator begin()  const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->__entries[this->n_entries - 1]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->__entries[this->n_entries - 1]); }

	// Add
	inline bool addData(const buffer_type& buffer){ return(this->addData(buffer.data(), buffer.size())); }
	bool addData(const char* const data, const U64 l_data){
		assert(l_data % sizeof(value_type) == 0);
		const size_t entries_adding = l_data / sizeof(value_type);

		// Check
		if(entries_adding + this->size() > this->capacity()){
			std::cerr << "could not fit!" << std::endl;
			return false;
		}

		U32 cumulative_position = 0;
		size_t start_position = this->size();
		for(size_t i = 0; i < entries_adding; ++i){
			new( &this->__entries[start_position + i] ) value_type( &data[cumulative_position] );
			cumulative_position += sizeof(value_type);
		}
		assert(cumulative_position == l_data);

		return true;
	}

	bool addEntry(const char* const data){
		// Check
		if(this->size() + 1 > this->capacity()){
			std::cerr << "could not fit!" << std::endl;
			return false;
		}

		new ( &this->__entries[this->size()] ) value_type( data );

		return true;
	}

protected:
	size_type  n_entries;
	size_type  n_capacity;
	pointer    __entries;
};

}

#endif /* TOMAHAWK_BASE_OUTPUT_CONTAINER_H_ */
