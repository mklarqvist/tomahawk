#ifndef BASIC_READER_H_
#define BASIC_READER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <cstddef>
#include <cstring>

#include "../support/MagicConstants.h" // for SILENT

namespace Tomahawk{

/*
 * Basic block-wise reader
 * Reads block_size bytes per iteration
 */
class reader {
	typedef char		  type;
	typedef type          value_type;
	typedef type         *pointer;
	typedef const type   *const_pointer;
	typedef type         &reference;
	typedef const type   &const_reference;
	typedef size_t        size_type;
	typedef ptrdiff_t     difference_type;

public:
	reader();
	reader(std::string input);
	reader(std::string input, const size_t block_size);
	virtual ~reader(){ delete[] this->buffer_; }
	virtual const_reference operator[](const size_t p) const{ return this->buffer_[p]; }
	virtual reference operator[](const size_t p){ return this->buffer_[p]; }

	class const_iterator{
		typedef const_iterator self_type;
		typedef type 		   value_type;
		typedef type		  &reference;
		typedef const type    &const_reference;
		typedef type 		  *pointer;
		typedef const type    *const_pointer;
		typedef int 	   	   difference_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(const_pointer ptr) : ptr_(ptr) { }
		virtual self_type operator++() { this->ptr_++; return *this; }
		virtual self_type operator++(int junk) { self_type i = *this; this->ptr_++; return i; }
		virtual const_reference operator*(){ return *ptr_; }
		virtual const_pointer operator->(){ return ptr_; }
		virtual bool operator==(const self_type& rhs) const{ return this->ptr_ == rhs.ptr_; }
		virtual bool operator!=(const self_type& rhs) const{ return this->ptr_ != rhs.ptr_; }
		virtual self_type operator=(const self_type& other) { this->ptr_ = other.ptr_; return *this; }

	private:
		const_pointer ptr_;
	};

	const_iterator begin() const{ return const_iterator(&(*this)[0]); }
	const_iterator end() const{ return const_iterator(&(*this)[this->end_]); }

	void clear(void){ this->end_ = 0; }
	bool empty(void) const{ return(this->end_ == 0); }
	size_t capacity(void) const{ return this->capacity_; }

	void capacity(uint32_t newSize){
		delete [] this->buffer_;
		this->capacity_ = newSize;
		this->buffer_ = new type[newSize];
	}

	void resize(void){
		type* old = this->buffer_;
		this->buffer_ = new type[this->capacity_ * 2];
		memcpy(this->buffer_, old, this->end_);
		this->capacity_ *= 2;
		delete [] old;
	}

	void resize(uint32_t newSize){
		uint32_t newSizeInternal = newSize;
		if(newSizeInternal < this->end_){
			this->end_ = newSizeInternal;
			return;
		}

		if(newSizeInternal < this->capacity_)
			return;

		type* old = this->buffer_;
		this->buffer_ = new type[newSize * 2];
		memcpy(this->buffer_, old, newSizeInternal);
		this->capacity_ = newSize * 2;
		delete [] old;
	}

	virtual const size_t size(void) const { return this->end_; } // Virtual -> allowed to overwrite in children classes
	bool good(void) const{ return this->stream_.good(); }

	virtual bool open(void);
	virtual bool open(std::string filename);
	void close(void);
	virtual bool read(void);
	virtual bool read(const uint32_t length);
	virtual bool readAppend(const uint32_t length);
	inline const uint64_t& filesize(void){ return this->filesize_; }
	inline uint64_t tellg(void){ return this->stream_.tellg(); }
	bool getLine(void); // Read until finding a new line into buffer
	bool getLine(std::string& data); // Read until finding a new line into string

public:
	std::string filename_;	// Input file name
	uint64_t filesize_;		// Input file size
	size_t block_size_;		// Size of block read each iteration
	size_t capacity_;		// Capacity of buffer
	size_t end_;			// End pointer of data in buffer
	std::ifstream stream_;	// Input stream
	pointer buffer_;		// Buffer
};

}

#endif /* BASIC_READER_H_ */
