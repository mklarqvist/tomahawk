#include "buffer.h"

namespace tomahawk {

twk_buffer_t::twk_buffer_t() :
	owns_data_(true),
	n_chars_(0),
	width_(0),
	iterator_position_(0),
	buffer_(nullptr)
{}

twk_buffer_t::twk_buffer_t(const uint64_t size) :
	owns_data_(true),
	n_chars_(0),
	width_(size),
	iterator_position_(0),
	buffer_(new char[size])
{}

twk_buffer_t::twk_buffer_t(char* target, const size_t length) :
	owns_data_(false),
	n_chars_(length),
	width_(length),
	iterator_position_(0),
	buffer_(target)
{}

twk_buffer_t::twk_buffer_t(const uint64_t size, char* target) :
	owns_data_(false),
	n_chars_(0),
	width_(size),
	iterator_position_(0),
	buffer_(target)
{}

twk_buffer_t::twk_buffer_t(const self_type& other) :
	owns_data_(other.owns_data_),
	n_chars_(other.n_chars_),
	width_(other.width_),
	iterator_position_(other.iterator_position_),
	buffer_(new char[other.width_])
{
	memcpy(this->buffer_, other.buffer_, other.size());
}

twk_buffer_t::~twk_buffer_t(){
	if(this->owns_data_)
		delete [] this->buffer_;
}

twk_buffer_t::twk_buffer_t(self_type&& other) noexcept :
	owns_data_(other.owns_data_),
	n_chars_(other.n_chars_),
	width_(other.width_),
	iterator_position_(other.iterator_position_),
	buffer_(nullptr)
{
	std::swap(buffer_, other.buffer_);
	other.reset();
	other.width_ = 0;
}

twk_buffer_t& twk_buffer_t::operator=(const self_type& other){
	if(this->owns_data_) delete [] this->buffer_;
	this->owns_data_  = other.owns_data_;
	this->n_chars_    = other.n_chars_;
	this->width_      = other.width_;
	this->iterator_position_ = other.iterator_position_;
	this->buffer_     = new char[other.width_];
	memcpy(this->buffer_, other.buffer_, other.size());
	return(*this);
}

twk_buffer_t& twk_buffer_t::operator=(self_type&& other) noexcept{
	if(this->owns_data_) delete [] this->buffer_;
	this->buffer_     = nullptr;
	this->owns_data_  = other.owns_data_;
	this->n_chars_    = other.n_chars_;
	this->width_      = other.width_;
	this->iterator_position_ = other.iterator_position_;
	std::swap(this->buffer_, other.buffer_);
	other.reset();
	other.width_   = 0;
	return(*this);
}

void twk_buffer_t::resize(const uint64_t new_size){
	if(this->n_chars_ == 0 && new_size == 0) return;
	if(new_size < this->capacity()){
		if(this->n_chars_ > new_size)
			this->n_chars_ = new_size;
		return;
	}

	char* temp = new char[new_size];
	assert(this->size() < new_size);
	memcpy(temp, this->buffer_, this->size());
	delete [] this->buffer_;
	this->buffer_ = temp;
	this->width_  = new_size;
}

void twk_buffer_t::resize(const self_type& other){
	if(other.size() >= this->capacity()){
		this->resize(other.capacity());
	}
}

void twk_buffer_t::set(const size_t size){
	this->n_chars_ = 0;
	this->width_ = size;
	if(this->buffer_ != nullptr)
		delete [] this->buffer_;

	this->buffer_ = new char[size];
}

void twk_buffer_t::set(const size_t size, char* target){
	this->n_chars_ = 0;
	this->width_ = size;
	this->buffer_ = target;
}

void twk_buffer_t::set(char* target){
	this->n_chars_ = 0;
	this->width_ = 0;
	this->buffer_ = target;
}

void twk_buffer_t::Add(const char* data, const uint32_t length){
	if(this->size() + length >= this->capacity())
		this->resize((this->size() + length) * 2);

	memcpy(&this->buffer_[this->n_chars_], &data[0], length);
	this->n_chars_ += length;
}

void twk_buffer_t::AddReadble(const int8_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%d", value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const int16_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%d", value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const int32_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%d", value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const int64_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%" PRId64, value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const uint8_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%u", value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const uint16_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%u", value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const uint32_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%u", value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const uint64_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%" PRIu64, value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const float& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%g", value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const double& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%g", value);
	this->n_chars_ += ret;
}

void twk_buffer_t::AddReadble(const std::string& value){
	if(this->n_chars_ + value.size() >= this->width_)
		this->resize(std::max(this->n_chars_ + value.size() + 100, this->width_*2));
	*this += value;
}

twk_buffer_t& twk_buffer_t::operator+=(const self_type& other){
	if(this->size() + other.size() >= this->capacity())
		this->resize((this->size() + other.size()) * 2);

	memcpy(&this->buffer_[this->n_chars_], other.buffer_, other.n_chars_);
	this->n_chars_ += other.n_chars_;

	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const char& value){
	if(this->n_chars_ + sizeof(char) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	this->buffer_[this->n_chars_] = value;
	++this->n_chars_;
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const int8_t& value){
	if(this->n_chars_ + sizeof(int8_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	this->buffer_[this->n_chars_] = value;
	++this->n_chars_;
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const uint8_t& value){
	if(this->n_chars_ + sizeof(uint8_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	uint8_t* p = reinterpret_cast<uint8_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(uint8_t);
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const float& value){
	if(this->n_chars_ + sizeof(float) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	float* p = reinterpret_cast<float*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(float);
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const uint16_t& value){
	if(this->n_chars_ + sizeof(uint16_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	uint16_t* p = reinterpret_cast<uint16_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(uint16_t);
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const int16_t& value){
	if(this->n_chars_ + sizeof(int16_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	short* p = reinterpret_cast<short*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(short);
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const uint32_t& value){
	if(this->n_chars_ + sizeof(uint32_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	uint32_t* p = reinterpret_cast<uint32_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(uint32_t);
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const int32_t& value){
	if(this->n_chars_ + sizeof(int32_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	int32_t* p = reinterpret_cast<int32_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(int32_t);
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const double& value){
	if(this->n_chars_ + sizeof(double) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	double* p = reinterpret_cast<double*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(double);
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const uint64_t& value){
	if(this->n_chars_ + sizeof(uint64_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	uint64_t* p = reinterpret_cast<uint64_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(uint64_t);
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const int64_t& value){
	if(this->n_chars_ + sizeof(int64_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	int64_t* p = reinterpret_cast<int64_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(int64_t);
	return *this;
}

twk_buffer_t& twk_buffer_t::operator+=(const std::string& value){
	if(this->n_chars_ + value.size() + sizeof(uint8_t) >= this->width_){
		uint64_t resize_to = std::max(this->n_chars_ + value.size() + sizeof(uint8_t) + 1000, this->width_ * 2);
		this->resize(resize_to);
	}

	for(uint32_t i = 0; i < value.size(); ++i){
		this->buffer_[this->n_chars_] = value[i];
		++this->n_chars_;
	}

	return *this;
}

twk_buffer_t& operator>>(twk_buffer_t& data, uint8_t& target){
	target = *reinterpret_cast<uint8_t*>(&data.buffer_[data.iterator_position_++]);
	return(data);
}

twk_buffer_t& operator>>(twk_buffer_t& data, uint16_t& target){
	target = *reinterpret_cast<uint16_t*>(&data.buffer_[data.iterator_position_]);
	data.iterator_position_ += sizeof(uint16_t);
	return(data);
}

twk_buffer_t& operator>>(twk_buffer_t& data, uint32_t& target){
	target = *reinterpret_cast<uint32_t*>(&data.buffer_[data.iterator_position_]);
	data.iterator_position_ += sizeof(uint32_t);
	return(data);
}

twk_buffer_t& operator>>(twk_buffer_t& data, uint64_t& target){
	target = *reinterpret_cast<uint64_t*>(&data.buffer_[data.iterator_position_]);
	data.iterator_position_ += sizeof(uint64_t);
	return(data);
}

twk_buffer_t& operator>>(twk_buffer_t& data, int8_t& target){
	target = *reinterpret_cast<int8_t*>(&data.buffer_[data.iterator_position_++]);
	return(data);
}

twk_buffer_t& operator>>(twk_buffer_t& data, int16_t& target){
	target = *reinterpret_cast<int16_t*>(&data.buffer_[data.iterator_position_]);
	data.iterator_position_ += sizeof(int16_t);
	return(data);
}

twk_buffer_t& operator>>(twk_buffer_t& data, int32_t& target){
	target = *reinterpret_cast<int32_t*>(&data.buffer_[data.iterator_position_]);
	data.iterator_position_ += sizeof(int32_t);
	return(data);
}

twk_buffer_t& operator>>(twk_buffer_t& data, int64_t& target){
	target = *reinterpret_cast<int64_t*>(&data.buffer_[data.iterator_position_]);
	data.iterator_position_ += sizeof(int64_t);
	return(data);
}

twk_buffer_t& operator>>(twk_buffer_t& data, float& target){
	target = *reinterpret_cast<float*>(&data.buffer_[data.iterator_position_]);
	data.iterator_position_ += sizeof(float);
	return(data);
}

twk_buffer_t& operator>>(twk_buffer_t& data, double& target){
	target = *reinterpret_cast<double*>(&data.buffer_[data.iterator_position_]);
	data.iterator_position_ += sizeof(double);
	return(data);
}

std::ostream& operator<<(std::ostream& out, const twk_buffer_t& data){
	out.write(data.data(), data.size());
	return(out);
}

void SerializeString(const std::string& string, twk_buffer_t& buffer){
	uint32_t size_helper = string.size();
	SerializePrimitive(size_helper, buffer);
	buffer += string;
}

void DeserializeString(std::string& string, twk_buffer_t& buffer){
	uint32_t size_helper = 0;
	DeserializePrimitive(size_helper, buffer);
	string.resize(size_helper);
	buffer.read(&string[0], size_helper);
}

}
