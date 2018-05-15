#include "index_container.h"

namespace tomahawk{
namespace totempole{

IndexContainer::IndexContainer(void) :
	n_entries_(0),
	n_capacity_(1000),
	entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{

}

IndexContainer::IndexContainer(const size_t n_capacity_) :
	n_entries_(0),
	n_capacity_(n_capacity_),
	entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{

}

// Functions for when interpreting from a byte stream
// first value is the number of indices
IndexContainer::IndexContainer(const char* const data_buffer, const U32 l_data) :
	n_entries_(*reinterpret_cast<const size_type* const>(data_buffer)),
	n_capacity_(this->n_entries_),
	entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{
	U32 cumulative_position = sizeof(size_type);
	for(U32 i = 0; i < this->size(); ++i){
		new( &this->entries_[i] ) value_type(&data_buffer[cumulative_position]);
		cumulative_position += TWK_INDEX_ENTRY_SIZE;
	}
	assert(cumulative_position == l_data);
}

IndexContainer::~IndexContainer(){
	for(size_type i = 0; i < this->size(); ++i)
		((this->entries_ + i)->~IndexEntry)();

	::operator delete[](static_cast<void*>(this->entries_));
}

IndexContainer& IndexContainer::operator+=(const value_type& index_entry){
	if(this->size() + 1 >= this->capacity()){
		//std::cerr << "is full resizing" << std::endl;
		this->resize();
	}

	//std::cerr << helpers::timestamp("DEBUG") << "Adding: " << this->size() << "/" << this->capacity() << std::endl;
	new( &this->entries_[this->n_entries_] ) value_type(index_entry); // invoke copy ctor
	++this->n_entries_;
	return(*this);
}

void IndexContainer::resize(const size_t new_capacity){
	//std::cerr << helpers::timestamp("DEBUG") << "Resize: " << this->capacity() << "->" << new_capacity << std::endl;
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

std::pair<U32, U32> IndexContainer::findOverlap(const S32& contigID) const{
	// Find first hit
	size_t i = 0;
	for(; i < this->size(); ++i){
		if(this->at(i).overlaps(contigID))
			break;
	}

	if(i == this->size())
		return(std::pair<U32, U32>(0, 0));

	const size_t from = i;
	for(; i < this->size(); ++i){
		if(this->at(i).overlaps(contigID) == false)
			break;
	}

	return(std::pair<U32,U32>(from, i));
}

std::pair<U32, U32> IndexContainer::findOverlap(const S32& contigID, const U64& position) const{
	// Find first hit
	size_t i = 0;
	for(; i < this->size(); ++i){
		if(this->at(i).overlaps(contigID, position))
			break;
	}

	if(i == this->size())
		return(std::pair<U32, U32>(0, 0));

	const size_t from = i;
	for(; i < this->size(); ++i){
		if(this->at(i).overlaps(contigID, position) == false)
			break;
	}

	return(std::pair<U32,U32>(from, i));
}

std::pair<U32, U32> IndexContainer::findOverlap(const S32& contigID, const U64& from_position, const U64& to_position) const{
	// Find first hit
	size_t i = 0;
	for(; i < this->size(); ++i){
		if(this->at(i).overlaps(contigID, from_position, to_position))
			break;
	}

	if(i == this->size())
		return(std::pair<U32, U32>(0, 0));

	const size_t from = i;
	for(; i < this->size(); ++i){
		if(this->at(i).overlaps(contigID, from_position, to_position) == false)
			break;
	}

	return(std::pair<U32,U32>(from, i));
}

}
}
