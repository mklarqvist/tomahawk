#ifndef BCFENTRY_CPP_
#define BCFENTRY_CPP_

#include "BCFEntry.h"
#include "../../support/MagicConstants.h"

namespace Tomahawk {
namespace BCF {

BCFEntry::BCFEntry(void):
	pointer(0),
	limit(262144),
	l_ID(0),
	p_genotypes(0),
	ref_alt(0),
	data(new char[this->limit]),
	body(reinterpret_cast<body_type*>(this->data)),
	alleles(new string_type[100]),
	ID(nullptr),
	genotypes(nullptr)
{

}

BCFEntry::~BCFEntry(void){ delete [] this->data; }

void BCFEntry::resize(const U32 size){
	char* temp = this->data;
	this->data = new char[size];
	memcpy(this->data, temp, this->pointer);
	std::swap(temp, this->data);
	delete [] temp;
	this->body = reinterpret_cast<body_type*>(this->data);

	if(size > this->limit)
		this->limit = size;
}

void BCFEntry::add(const char* const data, const U32 length){
	if(this->pointer + length > this-> capacity())
		this->resize(this->pointer + length + 65536);

	memcpy(&this->data[this->pointer], data, length);
	this->pointer += length;
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
		const base_type& length = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
#if BCF_ASSERT == 1
		assert(length.low == 7);
#endif

		S32 finalLength = 0;
		switch(length.low){
		case(1): finalLength = *reinterpret_cast<const SBYTE* const>(&this->data[internal_pos++]); break;
		case(2): finalLength = *reinterpret_cast<const S16* const>(&this->data[internal_pos+=2]); break;
		case(3): finalLength = *reinterpret_cast<const S32* const>(&this->data[internal_pos+=4]); break;
		}
		this->l_ID = finalLength;
		this->ID = &this->data[internal_pos];

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

		this->alleles[i].length = alelle_base.high;
		this->alleles[i].data = &this->data[internal_pos];

		if(alelle_base.high == 15){
			const base_type& length = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
#if BCF_ASSERT == 1
		assert(length.low == 7);
#endif

			S32 finalLength = 0;
			switch(length.low){
			case(1): finalLength = *reinterpret_cast<const SBYTE* const>(&this->data[internal_pos++]); break;
			case(2): finalLength = *reinterpret_cast<const S16* const>(&this->data[internal_pos+=2]); break;
			case(3): finalLength = *reinterpret_cast<const S32* const>(&this->data[internal_pos+=4]); break;
			}

			this->alleles[i].length = finalLength;
		}
		internal_pos += this->alleles[i].length;
	}
}

bool BCFEntry::parse(void){
	U32 internal_pos = sizeof(body_type);
	this->__parseID(internal_pos);
	this->__parseRefAlt(internal_pos);
	this->SetRefAlt();

	internal_pos = this->body->l_shared + sizeof(U32)*2;
	const base_type& fmt_key = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
	const SBYTE& fmt_key_value = *reinterpret_cast<SBYTE*>(&this->data[internal_pos]);

	switch(fmt_key.low){
	case(1): case(7): ++internal_pos; break;
	case(2): internal_pos += 2; break;
	case(3): case(5): internal_pos += 4; break;
	}

	const base_type& fmt_type = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
	//std::cerr << "fmt_key:" << (int)fmt_key_value << '\t' <<  "fmt_type: " << (int)fmt_type.high << '\t' << (int)fmt_type.low << std::endl;

	this->genotypes = &this->data[internal_pos];
	this->p_genotypes = internal_pos;

	return true;
}

void BCFEntry::SetRefAlt(void){
	this->ref_alt = 0;

	switch(this->alleles[0].data[0]){
	case 'A': this->ref_alt ^= Tomahawk::Constants::REF_ALT_A << 4; break;
	case 'T': this->ref_alt ^= Tomahawk::Constants::REF_ALT_T << 4; break;
	case 'G': this->ref_alt ^= Tomahawk::Constants::REF_ALT_G << 4; break;
	case 'C': this->ref_alt ^= Tomahawk::Constants::REF_ALT_C << 4; break;
	case '.': this->ref_alt ^= Tomahawk::Constants::REF_ALT_N << 4; break;
	}

	switch(this->alleles[1].data[0]){
	case 'A': this->ref_alt ^= Tomahawk::Constants::REF_ALT_A << 0; break;
	case 'T': this->ref_alt ^= Tomahawk::Constants::REF_ALT_T << 0; break;
	case 'G': this->ref_alt ^= Tomahawk::Constants::REF_ALT_G << 0; break;
	case 'C': this->ref_alt ^= Tomahawk::Constants::REF_ALT_C << 0; break;
	case '.': this->ref_alt ^= Tomahawk::Constants::REF_ALT_N << 0; break;
	}
}

}
}



#endif /* BCFENTRY_CPP_ */
