#ifndef BASICBUFFER_H_
#define BASICBUFFER_H_

#include <cstddef>
#include <iostream>
#include "../support/TypeDefinitions.h"
#include "../support/helpers.h"

namespace Tomahawk {
namespace IO{

struct BasicBuffer{
	typedef BasicBuffer self_type;

	BasicBuffer() : pointer(0), width(0), data(nullptr){}
	BasicBuffer(const U64 size) : pointer(0), width(size), data(new char[size]){}
	BasicBuffer(char* target, const size_t length) : pointer(length), width(length), data(target){}
	BasicBuffer(const U64 size, char* target) : pointer(0), width(size), data(target){}
	BasicBuffer(const self_type& other) : pointer(0), width(other.width), data(new char[other.width]){}
	virtual ~BasicBuffer(){}

	inline void set(const size_t size){
		this->pointer = 0;
		this->width = size;
		if(this->data != nullptr)
			delete [] this->data;

		this->data = new char[size];
	}

	inline void deleteAll(void){ delete [] this->data; } // manual cleaup

	inline void set(const size_t size, char* target){
		this->pointer = 0;
		this->width = size;
		this->data = target;
	}

	inline virtual void set(char* target){
		this->pointer = 0;
		this->width = 0;
		this->data = target;
	}

	inline void reset(){ this->pointer = 0; }
	inline void move(const U64 to){ this->pointer = to; }
	inline const U64& size(void) const{ return this->pointer; }
	inline const U64& capacity(void) const{ return this->width; }

	void resize(const U64 new_size){
		if(new_size <= this->capacity()){
			if(new_size < this->size())
				this->pointer = new_size;

			return;
		}

		U64 copy_to = this->size();
		if(new_size < this->size()){
			copy_to = new_size;
			this->pointer = copy_to;
		}

		//std::cerr << Helpers::timestamp("DEBUG") << "Resizing buffer: " << this->capacity() << " -> " << new_size << "\tcopyto: " << copy_to << std::endl;
		char* target = this->data;
		this->data = new char[new_size];
		memcpy(&this->data[0], &target[0], copy_to);
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

		memcpy(&this->data[this->pointer], &data[0], length);
		this->pointer += length;
	}

	inline self_type& operator+=(const self_type& other){
		if(this->size() + other.size() >= this->capacity())
			this->resize((this->size() + other.size()) * 2);

		memcpy(&this->data[this->pointer], other.data, other.pointer);
		this->pointer += other.pointer;

		return *this;
	}

	inline self_type& operator+=(const char& value){
		if(this->pointer + sizeof(char) >= this->width)
			this->resize(this->width*2);

		this->data[this->pointer] = value;
		++this->pointer;
		return *this;
	}

	inline self_type& operator+=(const BYTE& value){
		if(this->pointer + sizeof(BYTE) >= this->width)
			this->resize(this->width*2);

		BYTE* p = reinterpret_cast<BYTE*>(&this->data[this->pointer]);
		*p = value;
		this->pointer += sizeof(BYTE);
		return *this;
	}

	inline self_type& operator+=(const float& value){
		if(this->pointer + sizeof(float) >= this->width)
			this->resize(this->width*2);

		float* p = reinterpret_cast<float*>(&this->data[this->pointer]);
		*p = value;
		this->pointer += sizeof(float);
		return *this;
	}

	inline self_type& operator+=(const U16 value){
		if(this->pointer + sizeof(U16) >= this->width)
			this->resize(this->width*2);

		U16* p = reinterpret_cast<U16*>(&this->data[this->pointer]);
		*p = value;
		this->pointer += sizeof(U16);
		return *this;
	}

	inline self_type& operator+=(const short& value){
		if(this->pointer + sizeof(short) >= this->width)
			this->resize(this->width*2);

		short* p = reinterpret_cast<short*>(&this->data[this->pointer]);
		*p = value;
		this->pointer += sizeof(short);
		return *this;
	}

	inline self_type& operator+=(const U32& value){
		if(this->pointer + sizeof(U32) >= this->width)
			this->resize(this->width*2);

		U32* p = reinterpret_cast<U32*>(&this->data[this->pointer]);
		*p = value;
		this->pointer += sizeof(U32);
		return *this;
	}

	inline self_type& operator+=(const S32& value){
		if(this->pointer + sizeof(S32) >= this->width)
			this->resize(this->width*2);

		S32* p = reinterpret_cast<S32*>(&this->data[this->pointer]);
		*p = value;
		this->pointer += sizeof(S32);
		return *this;
	}

	inline self_type& operator+=(const double& value){
		if(this->pointer + sizeof(double) >= this->width)
			this->resize(this->width*2);

		double* p = reinterpret_cast<double*>(&this->data[this->pointer]);
		*p = value;
		this->pointer += sizeof(double);
		return *this;
	}

	inline self_type& operator+=(const U64& value){
		if(this->pointer + sizeof(U64) >= this->width)
			this->resize(this->width*2);

		U64* p = reinterpret_cast<U64*>(&this->data[this->pointer]);
		*p = value;
		this->pointer += sizeof(U64);
		return *this;
	}

	inline self_type& operator+=(const std::string& value){
		if(this->pointer + value.size() + sizeof(BYTE) >= this->width){
			U64 resize_to = this->width * 2;
			while(this->pointer + value.size() + sizeof(BYTE) >= resize_to)
				resize_to *= 2;

			this->resize(resize_to);
		}

		for(U32 i = 0; i < value.size(); ++i){
			this->data[this->pointer] = value[i];
			++this->pointer;
		}

		return *this;
	}


	char& operator[](const U64 size){ return this->data[size]; }
	const char& operator[](const U64 size) const{ return this->data[size]; }
	friend std::ostream& operator<<(std::ostream& out, const self_type& data){
		out.write(data.data, data.pointer);
		return(out);
	}

	U64 pointer;
	U64 width;
	char* data;
};

} /* namespace IO */
} /* namespace Tomahawk */

#endif /* BASICBUFFER_H_ */
