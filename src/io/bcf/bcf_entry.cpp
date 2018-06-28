#ifndef BCFENTRY_CPP_
#define BCFENTRY_CPP_

#include <support/magic_constants.h>
#include "bcf_entry.h"

#include <cassert>


namespace tomahawk {
namespace bcf {

BCFEntry::BCFEntry(void):
	l_data(0),
	l_capacity(262144),
	l_ID(0),
	ref_alt(0),
	isGood(false),
	data(new char[this->l_capacity]),
	body(reinterpret_cast<body_type*>(this->data)),
	alleles(new string_type[100]),
	ID(nullptr),
	hasGenotypes(false),
	ploidy(0),
	filter_start(0),
	n_filter(0),
	filterPointer(0),
	infoPointer(0),
	formatPointer(0),
	filterID(new BCFKeyTuple[256]),
	infoID(new BCFKeyTuple[256]),
	formatID(new BCFKeyTuple[256])
{

}

BCFEntry::BCFEntry(const self_type& other):
	l_data(other.l_data),
	l_capacity(other.l_capacity),
	l_ID(other.l_ID),
	ref_alt(other.ref_alt),
	isGood(other.isGood),
	data(new char[other.l_capacity]),
	body(nullptr),
	alleles(new string_type[100]),
	ID(nullptr),
	hasGenotypes(other.hasGenotypes),
	ploidy(other.ploidy),
	filter_start(other.filter_start),
	n_filter(other.n_filter),
	filterPointer(other.filterPointer),
	infoPointer(other.infoPointer),
	formatPointer(other.formatPointer),
	gt_support(other.gt_support),
	filterID(new BCFKeyTuple[256]),
	infoID(new BCFKeyTuple[256]),
	formatID(new BCFKeyTuple[256])
{
	memcpy(this->data, other.data, other.l_data);
	this->body = reinterpret_cast<body_type*>(this->data);

	for(U32 i = 0; i < 256; ++i){
		this->filterID[i] = other.filterID[i];
		this->formatID[i] = other.formatID[i];
		this->infoID[i]   = other.infoID[i];
	}

	U32 internal_pos = sizeof(body_type);
	this->__parseID(internal_pos);
	for(U32 i = 0; i < this->body->n_allele; ++i)
		this->alleles[i] = other.alleles[i];
}

BCFEntry::BCFEntry(self_type&& other) noexcept :
	l_data(other.l_data),
	l_capacity(other.l_capacity),
	l_ID(other.l_ID),
	ref_alt(other.ref_alt),
	isGood(other.isGood),
	data(other.data),
	body(other.body),
	alleles(other.alleles),
	ID(other.ID),
	hasGenotypes(other.hasGenotypes),
	ploidy(other.ploidy),
	filter_start(other.filter_start),
	n_filter(other.n_filter),
	filterPointer(other.filterPointer),
	infoPointer(other.infoPointer),
	formatPointer(other.formatPointer),
	gt_support(other.gt_support),
	filterID(other.filterID),
	infoID(other.infoID),
	formatID(other.formatID)
{
	other.data      = nullptr;
	other.body      = nullptr;
	other.alleles   = nullptr;
	other.ID        = nullptr;
	other.filterID  = nullptr;
	other.infoID    = nullptr;
	other.formatID  = nullptr;
}

BCFEntry& BCFEntry::operator=(const self_type& other){
	self_type tmp(other);         // re-use copy-constructor
	*this = std::move(tmp); // re-use move-assignment
	return *this;
}

BCFEntry& BCFEntry::operator=(self_type&& other) noexcept{
	if (this == &other)
		return *this;

	delete [] this->data;
	delete [] this->alleles;
	delete [] this->filterID;
	delete [] this->infoID;
	delete [] this->formatID;
	this->l_data        = other.l_data;
	this->l_capacity    = other.l_capacity;
	this->l_ID          = other.l_ID;
	this->ref_alt       = other.ref_alt;
	this->isGood        = other.isGood;
	this->data          = other.data;
	this->body          = other.body;
	this->alleles       = other.alleles;
	this->ID            = other.ID;
	this->hasGenotypes  = other.hasGenotypes;
	this->ploidy        = other.ploidy;
	this->filter_start  = other.filter_start;
	this->n_filter      = other.n_filter;
	this->filterPointer = other.filterPointer;
	this->infoPointer   = other.infoPointer;
	this->formatPointer = other.formatPointer;
	this->filterID      = other.filterID;
	this->infoID        = other.infoID;
	this->formatID      = other.formatID;
	this->gt_support    = other.gt_support;

	other.data     = nullptr;
	other.alleles  = nullptr;
	other.filterID = nullptr;
	other.infoID   = nullptr;
	other.formatID = nullptr;

	return *this;
}

BCFEntry::~BCFEntry(void){
	delete [] this->data;
	delete [] this->alleles;
	delete [] this->filterID;
	delete [] this->infoID;
	delete [] this->formatID;
}

void BCFEntry::resize(const U32 size){
	if(size == 0)
		return;

	char* temp = this->data;
	this->data = new char[size];
	memcpy(this->data, temp, this->l_data);
	delete [] temp;
	this->body = reinterpret_cast<body_type*>(this->data);

	if(size > this->l_capacity)
		this->l_capacity = size;
}

void BCFEntry::add(const char* const data, const U32 length){
	if(this->l_data + length > this-> capacity())
		this->resize(this->l_data + length + 65536);

	memcpy(&this->data[this->l_data], data, length);
	this->l_data += length;
}

void BCFEntry::__parseID(U32& internal_pos){
	// Parse ID
	const base_type& ID_base = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
#if BCF_ASSERT == 1
	assert(ID_base.low == 7);
#endif

	this->ID = &this->data[internal_pos];
	this->l_ID = ID_base.high;
	if(ID_base.high == 0){ // has no name
		this->l_ID = 0;
	} else if(ID_base.high == 15){
		// next byte is the length array
		// Type and length
		const base_type& array_base = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
		this->l_ID = this->getInteger(array_base.low, internal_pos);
		this->ID   = &this->data[internal_pos];
	}
	internal_pos += this->l_ID;
}

void BCFEntry::__parseRefAlt(U32& internal_pos){
	// Parse REF-ALT
	for(U32 i = 0; i < this->body->n_allele; ++i){
		const base_type& alelle_base = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
#if BCF_ASSERT == 1
		assert(alelle_base.low == 7);
#endif

		S32 length =  alelle_base.high;
		const char* ref_alt_data  = &this->data[internal_pos];

		if(alelle_base.high == 15){
			// Type and length
			const base_type& array_base = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
			length =  this->getInteger(array_base.low, internal_pos);
			ref_alt_data = &this->data[internal_pos];
		}
		this->alleles[i](ref_alt_data, length);

		//std::cerr << std::string(this->alleles[i].data, this->alleles[i].length) << std::endl;
		internal_pos += this->alleles[i].length;
	}
}

bool BCFEntry::nextFilter(S32& value, U32& position){
	if(this->filterPointer == this->n_filter)
		return false;

	value = this->getInteger(this->filter_key.low, position);
	this->filterID[this->filterPointer++].mapID = value;

	return true;
}

bool BCFEntry::nextInfo(S32& value, U32& length, BYTE& value_type, U32& position){
	if(this->infoPointer == this->body->n_info)
		return false;

	const base_type& info_key = *reinterpret_cast<const base_type* const>(&this->data[position++]);
	#if BCF_ASSERT == 1
	// The first object returns a single identifier
	// to a field. It should always be a single
	// value
	assert(info_key.high == 1);
	#endif

	// INFO identifier
	value = this->getInteger(info_key.low, position);
	this->infoID[this->infoPointer++].mapID = value;

	// Data for this identifier
	const base_type& info_value = *reinterpret_cast<const base_type* const>(&this->data[position++]);
	length = info_value.high;
	if(length == 15){
		const base_type& array_base = *reinterpret_cast<const base_type* const>(&this->data[position++]);
		length = this->getInteger(array_base.low, position);
	}
	value_type = info_value.low;

	return true;
}

bool BCFEntry::nextFormat(S32& value, U32& length, BYTE& value_type, U32& position){
	if(this->formatPointer == this->body->n_fmt)
		return false;

	const base_type& format_key = *reinterpret_cast<const base_type* const>(&this->data[position++]);
	#if BCF_ASSERT == 1
	// This first bit returns a single identifier
	// to a field. It should always be a single
	// value
	assert(format_key.high == 1);
	#endif

	// format identifier
	value = this->getInteger(format_key.low, position);
	this->formatID[this->formatPointer++].mapID = value;

	// Data for this identifier
	const base_type& format_value = *reinterpret_cast<const base_type* const>(&this->data[position++]);
	length = format_value.high;
	if(length == 15){
		const base_type& array_base = *reinterpret_cast<const base_type* const>(&this->data[position++]);
		length = this->getInteger(array_base.low, position);
	}
	value_type = format_value.low;

	return true;
}

bool BCFEntry::parse(const U64 n_samples){
	this->body = reinterpret_cast<body_type*>(this->data);
	U32 internal_pos = sizeof(body_type);
	this->__parseID(internal_pos);
	this->__parseRefAlt(internal_pos);
	this->SetRefAlt();

	// start of FILTER
	this->filter_start = internal_pos;

	// At FILTER
	// Typed vector
	const base_type& filter_key = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
	U32 n_filter     = filter_key.high;
	if(n_filter == 15) n_filter = this->getInteger(filter_key.low, internal_pos);
	this->n_filter   = n_filter;
	this->filter_key = filter_key;

	S32 val = 0;
	while(this->nextFilter(val, internal_pos)){}

	// At INFO
	U32  info_length;
	BYTE info_primitive_type;
	for(U32 i = 0; i < this->body->n_info; ++i){
		if(this->nextInfo(val, info_length, info_primitive_type, internal_pos) == false){
			std::cerr << "illegal match info" << std::endl;
			exit(1);
		}
		this->infoID[i].l_stride       = info_length;
		this->infoID[i].primitive_type = info_primitive_type;
		this->infoID[i].l_offset       = internal_pos;

		// Flags and integers
		// These are BCF value types
		if(info_primitive_type <= 3){
			for(U32 j = 0; j < info_length; ++j){
				this->getInteger(info_primitive_type, internal_pos);
			}
		}
		// Floats
		else if(info_primitive_type == 5){
			for(U32 j = 0; j < info_length; ++j){
				this->getFloat(internal_pos);
			}
		}
		// Chars
		else if(info_primitive_type == 7){
			for(U32 j = 0; j < info_length; ++j){
				this->getChar(internal_pos);
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible in info: " << (int)info_primitive_type << std::endl;
			exit(1);
		}
	}

	assert(internal_pos == (this->body->l_shared + sizeof(U32)*2));

	if(internal_pos == this->l_data){
		//std::cerr << "have no FORMAT data" << std::endl;
		return true;
	}

	BYTE format_primitive_type = 0;
	for(U32 i = 0; i < this->body->n_fmt; ++i){
		if(this->nextFormat(val, info_length, format_primitive_type, internal_pos) == false){
			std::cerr << "illegal match format" << std::endl;
			exit(1);
		}

		this->formatID[i].l_stride       = info_length;
		this->formatID[i].primitive_type = format_primitive_type;
		this->formatID[i].l_offset       = internal_pos;

		// Flags and integers
		// These are BCF value types
		if(format_primitive_type <= 3){
			for(U32 s = 0; s < n_samples; ++s){
				for(U32 j = 0; j < info_length; ++j)
					this->getInteger(format_primitive_type, internal_pos);
			}
		}
		// Floats
		else if(format_primitive_type == 5){
			for(U32 s = 0; s < n_samples; ++s){
				for(U32 j = 0; j < info_length; ++j)
					this->getFloat(internal_pos);

			}
		}
		// Chars
		else if(format_primitive_type == 7){
			for(U32 s = 0; s < n_samples; ++s){
				for(U32 j = 0; j < info_length; ++j)
					this->getChar(internal_pos);
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible in format: " << (int)format_primitive_type << std::endl;
			std::cerr << helpers::timestamp("LOG") << val << '\t' << info_length << '\t' << (int)format_primitive_type << '\t' << internal_pos << '/' << this->l_data << std::endl;
			exit(1);
		}
	}

	assert(internal_pos == this->l_data);

	this->isGood = true;
	return true;
}

void BCFEntry::SetRefAlt(void){
	this->ref_alt = 0;

	switch(this->alleles[0].data[0]){
	case 'A': this->ref_alt ^= tomahawk::constants::REF_ALT_A << 4; break;
	case 'T': this->ref_alt ^= tomahawk::constants::REF_ALT_T << 4; break;
	case 'G': this->ref_alt ^= tomahawk::constants::REF_ALT_G << 4; break;
	case 'C': this->ref_alt ^= tomahawk::constants::REF_ALT_C << 4; break;
	case '.': this->ref_alt ^= tomahawk::constants::REF_ALT_N << 4; break;
	}

	switch(this->alleles[1].data[0]){
	case 'A': this->ref_alt ^= tomahawk::constants::REF_ALT_A << 0; break;
	case 'T': this->ref_alt ^= tomahawk::constants::REF_ALT_T << 0; break;
	case 'G': this->ref_alt ^= tomahawk::constants::REF_ALT_G << 0; break;
	case 'C': this->ref_alt ^= tomahawk::constants::REF_ALT_C << 0; break;
	case '.': this->ref_alt ^= tomahawk::constants::REF_ALT_N << 0; break;
	}
}

}
}



#endif /* BCFENTRY_CPP_ */
