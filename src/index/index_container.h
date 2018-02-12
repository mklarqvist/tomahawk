#ifndef INDEX_INDEX_CONTAINER_H_
#define INDEX_INDEX_CONTAINER_H_

#include <cstring>  // size_t, ptrdiff_t
#include <iterator> // forward_iterator_tag

#include "../support/type_definitions.h"
#include "../io/BasicBuffer.h"
#include "index_entry.h"

namespace Tomahawk{
namespace Totempole{

/**<
 * STL-like container for Tomahawk index entries
 */
class IndexContainer{
private:
	typedef IndexContainer        self_type;
	typedef IndexEntry            value_type;
    typedef value_type&           reference;
    typedef const value_type&     const_reference;
    typedef value_type*           pointer;
    typedef const value_type*     const_pointer;
    typedef std::ptrdiff_t        difference_type;
    typedef std::size_t           size_type;
    typedef IO::BasicBuffer       buffer_type;

public:
    IndexContainer(void) :
		n_entries_(0),
		n_capacity_(1000),
		entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
	{

	}

    IndexContainer(const size_t n_capacity_) :
    	n_entries_(0),
		n_capacity_(n_capacity_),
		entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
	{

	}

    // Functions for when interpreting from a byte stream
    // first value is the number of indices
	IndexContainer(const char* const data_buffer);
	IndexContainer(const buffer_type& data_buffer);

	~IndexContainer(){
		for(size_type i = 0; i < this->size(); ++i)
			((this->entries_ + i)->~IndexEntry)();

		::operator delete[](static_cast<void*>(this->entries_));
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
	inline reference at(const size_type& position){ return(this->entries_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->entries_[position]); }
	inline reference operator[](const size_type& position){ return(this->entries_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->entries_[position]); }
	inline pointer data(void){ return(this->entries_); }
	inline const_pointer data(void) const{ return(this->entries_); }
	inline reference front(void){ return(this->entries_[0]); }
	inline const_reference front(void) const{ return(this->entries_[0]); }
	inline reference back(void){ return(this->entries_[this->n_entries_ - 1]); }
	inline const_reference back(void) const{ return(this->entries_[this->n_entries_ - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries_ == 0); }
	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->entries_[0]); }
	inline iterator end()  { return iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator begin()  const{ return const_iterator(&this->entries_[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->entries_[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->entries_[this->n_entries_]); }

	// Overload basic operator
	self_type& operator+=(const value_type& index_entry){
		if(this->size() + 1 >= this->capacity()){
			//std::cerr << "is full resizing" << std::endl;
			this->resize();
		}

		//std::cerr << Helpers::timestamp("DEBUG") << "Adding: " << this->size() << "/" << this->capacity() << std::endl;
		new( &this->entries_[this->n_entries_] ) value_type(index_entry); // invoke copy ctor
		++this->n_entries_;
		return(*this);
	}

	void resize(const size_t new_capacity){
		//std::cerr << Helpers::timestamp("DEBUG") << "Resize: " << this->capacity() << "->" << new_capacity << std::endl;
		// if resizing to a smaller size
		if(new_capacity < this->capacity()){
			// Call destructor for values between shrunk size and previous numbers
			for(size_type i = new_capacity; i < this->size(); ++i)
				((this->entries_ + i)->~IndexEntry)();

			this->n_entries_ = new_capacity;
			return;
		}

		pointer temp = this->entries_; // Move current data pointer
		this->entries_ = static_cast<pointer>(::operator new[](new_capacity*sizeof(value_type))); // Allocate new memory at old pointer
		// Copy data over from temporary data pointer to new pointer
		for(U32 i = 0; i < this->size(); ++i)
			new( &this->entries_[i] ) value_type(temp[i]);

		// Release memory from the temporary address
		for(size_type i = 0; i < this->size(); ++i)
			((temp + i)->~IndexEntry)();

		::operator delete[](static_cast<void*>(temp));
		this->n_capacity_ = new_capacity;
	}
	inline void resize(void){ this->resize(this->capacity()*2); }

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& container){
		stream.write(reinterpret_cast<const char*>(&container.n_entries_), sizeof(size_type));
		for(size_type i = 0; i < container.size(); ++i)
			stream << container[i];

		return stream;
	}

private:
	size_type  n_entries_;
	size_type  n_capacity_;
	pointer    entries_;
};

}
}

#endif /* INDEX_INDEX_CONTAINER_H_ */
