/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TOMAHAWK_BUFFER_H_
#define TOMAHAWK_BUFFER_H_

#define __STDC_FORMAT_MACROS 1
#include <inttypes.h>

#include <cstddef>
#include <iostream>
#include <stdint.h>
#include <cassert>

#include "utility.h"
#include "generic_iterator.h"

namespace tomahawk {

struct twk_buffer_t {
public:
    typedef twk_buffer_t      self_type;
    typedef char              value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;
    typedef std::size_t       size_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
	twk_buffer_t();
	twk_buffer_t(const uint64_t size);
	twk_buffer_t(char* target, const size_t length);
	twk_buffer_t(const uint64_t size, char* target);
	twk_buffer_t(const self_type& other);
	twk_buffer_t(self_type&& other) noexcept;
	self_type& operator=(const self_type& other);
	self_type& operator=(self_type&& other) noexcept;
	~twk_buffer_t();

	inline reference back(void){ return(this->buffer_[this->n_chars_-1]); }
	inline reference front(void){ return(this->buffer_[0]); }
	inline const_reference back(void) const{ return(this->buffer_[this->n_chars_-1]); }
	inline const_reference front(void) const{ return(this->buffer_[0]); }

	// Iterator
	inline iterator begin(){ return iterator(&this->buffer_[0]); }
	inline iterator end()  { return iterator(&this->buffer_[this->n_chars_]); }
	inline const_iterator begin()  const{ return const_iterator(&this->buffer_[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->buffer_[this->n_chars_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->buffer_[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->buffer_[this->n_chars_]); }

	inline void set(const size_t size);
	inline void set(const size_t size, char* target);
	inline void set(char* target);

	inline void clear(){ this->n_chars_ = 0; this->width_ = 0; this->iterator_position_ = 0; if(owns_data_){ delete[] buffer_; buffer_ = nullptr; } }
	inline void reset(){ this->n_chars_ = 0; this->iterator_position_ = 0; }
	inline void resetIterator(){ this->iterator_position_ = 0; }
	inline void move(const uint64_t to){ this->n_chars_ = to; }
	inline const uint64_t& size(void) const{ return this->n_chars_; }
	inline const uint64_t& capacity(void) const{ return this->width_; }

	void resize(const uint64_t new_size);
	void resize(const self_type& other);

	void Add(const char* data, const uint32_t length);
	void AddReadble(const int8_t& value);
	void AddReadble(const int16_t& value);
	void AddReadble(const int32_t& value);
	void AddReadble(const int64_t& value);
	void AddReadble(const uint8_t& value);
	void AddReadble(const uint16_t& value);
	void AddReadble(const uint32_t& value);
	void AddReadble(const uint64_t& value);
	void AddReadble(const float& value);
	void AddReadble(const double& value);
	void AddReadble(const std::string& value);
	self_type& operator+=(const self_type& other);
	self_type& operator+=(const char& value);
	self_type& operator+=(const int8_t& value);
	self_type& operator+=(const uint8_t& value);
	self_type& operator+=(const float& value);
	self_type& operator+=(const uint16_t& value);
	self_type& operator+=(const int16_t& value);
	self_type& operator+=(const uint32_t& value);
	self_type& operator+=(const int32_t& value);
	self_type& operator+=(const double& value);
	self_type& operator+=(const uint64_t& value);
	self_type& operator+=(const int64_t& value);
	self_type& operator+=(const std::string& value);

	inline reference operator[](const uint64_t position){ return this->buffer_[position]; }
	inline const_reference operator[](const uint64_t position) const{ return this->buffer_[position]; }
	inline reference at(const uint64_t position){ return this->buffer_[position]; }
	inline const_reference at(const uint64_t position) const{ return this->buffer_[position]; }
	inline pointer data(void){ return(this->buffer_); }
	inline const_pointer data(void) const{ return(this->buffer_); }

	void read(char* target, const uint32_t n_length){
		memcpy(target, &this->buffer_[this->iterator_position_], n_length);
		this->iterator_position_ += n_length;
	}

private:
	friend self_type& operator>>(self_type& data, uint8_t& target);
	friend self_type& operator>>(self_type& data, uint16_t& target);
	friend self_type& operator>>(self_type& data, uint32_t& target);
	friend self_type& operator>>(self_type& data, uint64_t& target);
	friend self_type& operator>>(self_type& data, int8_t& target);
	friend self_type& operator>>(self_type& data, int16_t& target);
	friend self_type& operator>>(self_type& data, int32_t& target);
	friend self_type& operator>>(self_type& data, int64_t& target);
	friend self_type& operator>>(self_type& data, float& target);
	friend self_type& operator>>(self_type& data, double& target);

	friend std::ostream& operator<<(std::ostream& out, const self_type& data);

public:
	bool     owns_data_;
	uint64_t n_chars_;
	uint64_t width_;
	uint64_t iterator_position_;
	pointer  buffer_;
};

/**<
 * Supportive functions for serializing/deserialize data to/from a byte
 * stream.
 * @param value  Src value.
 * @param buffer Dst buffer reference.
 */
void SerializeString(const std::string& string, twk_buffer_t& buffer);
void DeserializeString(std::string& string, twk_buffer_t& buffer);

template <class T>
static inline void SerializePrimitive(const T& value, twk_buffer_t& buffer){
	buffer += value;
}

template <class T>
static inline void DeserializePrimitive(T& value, twk_buffer_t& buffer){
	buffer >> value;
}

}

#endif
