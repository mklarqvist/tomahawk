#ifndef BASICBUFFER_H_
#define BASICBUFFER_H_

#include <cstddef>
#include <iostream>
#include "support/type_definitions.h"
#include "support/helpers.h"

namespace tomahawk{
namespace io{

struct BasicBuffer{
private:
    typedef BasicBuffer       self_type;
    typedef char              value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;
    typedef std::size_t       size_type;

public:
	BasicBuffer() : owns_data_(true), n_chars(0), width(0), iterator_position_(0), buffer(nullptr){}
	BasicBuffer(const U64 size) : owns_data_(true), n_chars(0), width(size), iterator_position_(0), buffer(new char[size]){}
	BasicBuffer(char* target, const size_t length) : owns_data_(false), n_chars(length), width(length), iterator_position_(0), buffer(target){}
	BasicBuffer(const U64 size, char* target) : owns_data_(false), n_chars(0), width(size), iterator_position_(0), buffer(target){}
	BasicBuffer(const self_type& other) : owns_data_(other.owns_data_), n_chars(0), width(other.width), iterator_position_(0), buffer(new char[other.width]){}
	virtual ~BasicBuffer(){
		if(this->owns_data_)
			delete [] this->buffer;
	}

	inline reference back(void){ return(this->buffer[this->n_chars-1]); }
	inline reference front(void){ return(this->buffer[0]); }
	inline const_reference back(void) const{ return(this->buffer[this->n_chars-1]); }
	inline const_reference front(void) const{ return(this->buffer[0]); }

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

	// Iterator
	inline iterator begin(){ return iterator(&this->buffer[0]); }
	inline iterator end()  { return iterator(&this->buffer[this->n_chars]); }
	inline const_iterator begin()  const{ return const_iterator(&this->buffer[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->buffer[this->n_chars]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->buffer[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->buffer[this->n_chars]); }

	inline void set(const size_t size){
		this->n_chars = 0;
		this->width = size;
		if(this->buffer != nullptr)
			delete [] this->buffer;

		this->buffer = new char[size];
	}

	inline void set(const size_t size, char* target){
		this->n_chars = 0;
		this->width = size;
		this->buffer = target;
	}

	inline virtual void set(char* target){
		this->n_chars = 0;
		this->width = 0;
		this->buffer = target;
	}

	inline void reset(){ this->n_chars = 0; this->iterator_position_ = 0; }
	inline void resetIterator(){ this->iterator_position_ = 0; }
	inline void move(const U64 to){ this->n_chars = to; }
	inline const U64& size(void) const{ return this->n_chars; }
	inline const U64& capacity(void) const{ return this->width; }

	void resize(const U64 new_size){
		if(new_size <= this->capacity()){
			if(new_size < this->size())
				this->n_chars = new_size;

			return;
		}

		U64 copy_to = this->size();
		if(new_size < this->size()){
			copy_to = new_size;
			this->n_chars = copy_to;
		}

		//std::cerr << utility::timestamp("DEBUG") << "Resizing buffer: " << this->capacity() << " -> " << new_size << "\tcopyto: " << copy_to << std::endl;
		char* target = this->buffer;
		this->buffer = new char[new_size];
		memcpy(&this->buffer[0], &target[0], copy_to);
		delete [] target;
		this->width = new_size;
	}

	void resize(const self_type& other){
		if(other.size() >= this->capacity()){
			this->resize(other.capacity());
		}
	}

	void Add(const char* data, const U32 length){
		if(this->size() + length >= this->capacity())
			this->resize((this->size() + length) * 2);

		memcpy(&this->buffer[this->n_chars], &data[0], length);
		this->n_chars += length;
	}

	void AddReadble(const SBYTE& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(this->width*2);
		const int ret = sprintf(&this->buffer[this->n_chars], "%d", value);
		this->n_chars += ret;
	}

	void AddReadble(const S16& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(this->width*2);
		const int ret = sprintf(&this->buffer[this->n_chars], "%d", value);
		this->n_chars += ret;
	}

	void AddReadble(const S32& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(this->width*2);
		const int ret = sprintf(&this->buffer[this->n_chars], "%d", value);
		this->n_chars += ret;
	}

	void AddReadble(const BYTE& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(this->width*2);
		const int ret = sprintf(&this->buffer[this->n_chars], "%u", value);
		this->n_chars += ret;
	}

	void AddReadble(const U16& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(this->width*2);
		const int ret = sprintf(&this->buffer[this->n_chars], "%u", value);
		this->n_chars += ret;
	}

	void AddReadble(const U32& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(this->width*2);
		const int ret = sprintf(&this->buffer[this->n_chars], "%u", value);
		this->n_chars += ret;
	}

	void AddReadble(const U64& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(this->width*2);
		const int ret = sprintf(&this->buffer[this->n_chars], "%llu", value);
		this->n_chars += ret;
	}

	void AddReadble(const float& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(this->width*2);
		const int ret = sprintf(&this->buffer[this->n_chars], "%g", value);
		this->n_chars += ret;
	}

	void AddReadble(const double& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(this->width*2);
		const int ret = sprintf(&this->buffer[this->n_chars], "%g", value);
		this->n_chars += ret;
	}

	void AddReadble(const std::string& value){
		if(this->n_chars + value.size() >= this->width)
			this->resize(this->width*2);
		*this += value;
	}

	inline self_type& operator+=(const self_type& other){
		if(this->size() + other.size() >= this->capacity())
			this->resize((this->size() + other.size()) * 2);

		memcpy(&this->buffer[this->n_chars], other.buffer, other.n_chars);
		this->n_chars += other.n_chars;

		return *this;
	}

	inline self_type& operator+=(const char& value){
		if(this->n_chars + sizeof(char) >= this->width)
			this->resize(this->width*2);

		this->buffer[this->n_chars] = value;
		++this->n_chars;
		return *this;
	}

	inline self_type& operator+=(const BYTE& value){
		if(this->n_chars + sizeof(BYTE) >= this->width)
			this->resize(this->width*2);

		BYTE* p = reinterpret_cast<BYTE*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(BYTE);
		return *this;
	}

	inline self_type& operator+=(const float& value){
		if(this->n_chars + sizeof(float) >= this->width)
			this->resize(this->width*2);

		float* p = reinterpret_cast<float*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(float);
		return *this;
	}

	inline self_type& operator+=(const U16& value){
		if(this->n_chars + sizeof(U16) >= this->width)
			this->resize(this->width*2);

		U16* p = reinterpret_cast<U16*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(U16);
		return *this;
	}

	inline self_type& operator+=(const short& value){
		if(this->n_chars + sizeof(short) >= this->width)
			this->resize(this->width*2);

		short* p = reinterpret_cast<short*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(short);
		return *this;
	}

	inline self_type& operator+=(const U32& value){
		if(this->n_chars + sizeof(U32) >= this->width)
			this->resize(this->width*2);

		U32* p = reinterpret_cast<U32*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(U32);
		return *this;
	}

	inline self_type& operator+=(const S32& value){
		if(this->n_chars + sizeof(S32) >= this->width)
			this->resize(this->width*2);

		S32* p = reinterpret_cast<S32*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(S32);
		return *this;
	}

	inline self_type& operator+=(const double& value){
		if(this->n_chars + sizeof(double) >= this->width)
			this->resize(this->width*2);

		double* p = reinterpret_cast<double*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(double);
		return *this;
	}

	inline self_type& operator+=(const U64& value){
		if(this->n_chars + sizeof(U64) >= this->width)
			this->resize(this->width*2);

		U64* p = reinterpret_cast<U64*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(U64);
		return *this;
	}

	inline self_type& operator+=(const size_t& value){
		if(this->n_chars + sizeof(size_t) >= this->width)
			this->resize(this->width*2);

		size_t* p = reinterpret_cast<size_t*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(size_t);
		return *this;
	}

	inline self_type& operator+=(const std::string& value){
		if(this->n_chars + value.size() + sizeof(BYTE) >= this->width){
			U64 resize_to = this->width * 2;
			while(this->n_chars + value.size() + sizeof(BYTE) >= resize_to)
				resize_to *= 2;

			this->resize(resize_to);
		}

		for(U32 i = 0; i < value.size(); ++i){
			this->buffer[this->n_chars] = value[i];
			++this->n_chars;
		}

		return *this;
	}

	inline reference operator[](const U64 position){ return this->buffer[position]; }
	inline const_reference operator[](const U64 position) const{ return this->buffer[position]; }
	inline reference at(const U64 position){ return this->buffer[position]; }
	inline const_reference at(const U64 position) const{ return this->buffer[position]; }
	inline pointer data(void){ return(this->buffer); }
	inline const_pointer data(void) const{ return(this->buffer); }

	void read(char* target, const U32 n_length){
		memcpy(target, &this->buffer[this->iterator_position_], n_length);
		this->iterator_position_ += n_length;
	}

private:
	friend self_type& operator>>(self_type& data, BYTE& target){
		target = *reinterpret_cast<BYTE*>(&data.buffer[data.iterator_position_++]);
		return(data);
	}

	friend self_type& operator>>(self_type& data, U16& target){
		target = *reinterpret_cast<U16*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(U16);
		return(data);
	}

	friend self_type& operator>>(self_type& data, U32& target){
		target = *reinterpret_cast<U32*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(U32);
		return(data);
	}

	friend self_type& operator>>(self_type& data, U64& target){
		target = *reinterpret_cast<U64*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(U64);
		return(data);
	}

	friend self_type& operator>>(self_type& data, SBYTE& target){
		target = *reinterpret_cast<SBYTE*>(&data.buffer[data.iterator_position_++]);
		return(data);
	}

	friend self_type& operator>>(self_type& data, S16& target){
		target = *reinterpret_cast<S16*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(S16);
		return(data);
	}

	friend self_type& operator>>(self_type& data, S32& target){
		target = *reinterpret_cast<S32*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(S32);
		return(data);
	}

	friend self_type& operator>>(self_type& data, int64_t& target){
		target = *reinterpret_cast<int64_t*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(int64_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, float& target){
		target = *reinterpret_cast<float*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(float);
		return(data);
	}

	friend self_type& operator>>(self_type& data, double& target){
		target = *reinterpret_cast<double*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(double);
		return(data);
	}

	friend self_type& operator>>(self_type& data, size_type& target){
		target = *reinterpret_cast<size_type*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(size_type);
		return(data);
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& data){
		out.write(data.data(), data.size());
		return(out);
	}

public:
	bool    owns_data_;
	U64     n_chars;
	U64     width;
	U64     iterator_position_;
	pointer buffer;
};

} /* namespace IO */
} /* namespace Tomahawk */

#endif /* BASICBUFFER_H_ */
