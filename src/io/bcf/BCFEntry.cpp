#ifndef BCFENTRY_CPP_
#define BCFENTRY_CPP_

#include <cassert>

#include "BCFEntry.h"
#include "../../support/MagicConstants.h"

namespace Tomahawk {
namespace BCF {

BCFEntry::BCFEntry(void):
	l_data(0),
	limit(262144),
	l_ID(0),
	p_genotypes(0),
	ref_alt(0),
	isGood(false),
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
	memcpy(this->data, temp, this->l_data);
	std::swap(temp, this->data);
	delete [] temp;
	this->body = reinterpret_cast<body_type*>(this->data);

	if(size > this->limit)
		this->limit = size;
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
		const base_type& length = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
#if BCF_ASSERT == 1
		assert(length.low == 7);
#endif

		S32 finalLength = 0;
		switch(length.low){
		case(1): finalLength = *reinterpret_cast<const SBYTE* const>(&this->data[internal_pos++]); break;
		case(2): finalLength = *reinterpret_cast<const S16* const>(&this->data[internal_pos+=2]);  break;
		case(3): finalLength = *reinterpret_cast<const S32* const>(&this->data[internal_pos+=4]);  break;
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
			case(2): finalLength = *reinterpret_cast<const S16* const>(&this->data[internal_pos+=2]);  break;
			case(3): finalLength = *reinterpret_cast<const S32* const>(&this->data[internal_pos+=4]);  break;
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

	// Format key
	const base_type& fmt_type = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
	//std::cerr << "fmt_key:" << (int)fmt_key_value << '\t' <<  "fmt_type: " << (int)fmt_type.high << '\t' << (int)fmt_type.low << std::endl;
	//std::cerr << (int)fmt_type_value2 << '\t' << (int)fmt_type_value1 << std::endl;
	//assert(fmt_type.high == 2);

	if(fmt_type.high != 2){
		this->isGood = false;
		return false;
	}

	this->isGood = true;

	/*
	for(U32 i = 0; i < 44; ++i){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
		std::cerr << i << ':' << " " << (int)fmt_type_value1 << ',' << (int)fmt_type_value2 << '\t' << (int)(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1)) << ',' << (int)(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2)) << std::endl;
	}
	*/
	this->genotypes = &this->data[internal_pos];
	this->p_genotypes = internal_pos;

	return true;
}

double BCFEntry::getMissingness(const U64& samples) const{
	if(!this->good())
		return(1);

	U32 internal_pos = this->p_genotypes;
	U64 n_missing = 0;
	for(U32 i = 0; i < samples; ++i){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
		//std::cerr << i << ':' << " " << (int)fmt_type_value1 << ',' << (int)fmt_type_value2 << '\t' << (int)(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1)) << ',' << (int)(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2)) << std::endl;

		if(fmt_type_value1 < 0 || fmt_type_value2 < 0)
			return(2);

		if(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) == 2 || BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) == 2) ++n_missing;
	}
	return((double)n_missing/samples);
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
