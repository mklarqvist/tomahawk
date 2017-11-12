#include <cassert>

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
	isGood(false),
	data(new char[this->limit]),
	body(reinterpret_cast<body_type*>(this->data)),
	alleles(new string_type[100]),
	ID(nullptr),
	genotypes(nullptr),
	filter_start(0),
	n_filter(0),
	// Vectors of identifiers
	filterPointer(0),
	infoPointer(0),
	formatPointer(0),
	// FILTER
	filterID(new U32[256]),
	// INFO
	infoID(new U32[256]),
	// FORMAT
	formatID(new U32[256])
{

}

BCFEntry::~BCFEntry(void){
	delete [] this->data;
	delete [] this->filterID;
	delete [] this->infoID;
	delete [] this->formatID;
	delete [] this->alleles;
}

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
		// Type and length
		const base_type& array_base = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
		this->l_ID = this->getInteger(array_base.low, internal_pos);
		this->ID = &this->data[internal_pos];
	}
	//std::cerr << std::string(this->ID, this->l_ID) << ":" << this->l_ID << '\t';
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
			// Type and length
			const base_type& array_base = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
			this->alleles[i].length = this->getInteger(array_base.low, internal_pos);
			this->alleles[i].data = &this->data[internal_pos];
		}

		//std::cerr << std::string(this->alleles[i].data, this->alleles[i].length) << '\t';
		internal_pos += this->alleles[i].length;
	}
}

bool BCFEntry::nextFilter(S32& value, U32& position){
	if(this->filterPointer == this->n_filter)
		return false;

	value = this->getInteger(this->filter_key.low, position);
	this->filterID[this->filterPointer++] = value;

	return true;
}

bool BCFEntry::nextInfo(S32& value, U32& length, BYTE& value_type, U32& position){
	if(this->infoPointer == this->body->n_info)
		return false;

	const base_type& info_key = *reinterpret_cast<const base_type* const>(&this->data[position++]);
	#if BCF_ASSERT == 1
	// This first bit returns a single identifier
	// to a field. It should always be a single
	// value
	assert(info_key.high == 1);
	#endif

	// INFO identifier
	value = this->getInteger(info_key.low, position);
	this->infoID[this->infoPointer++] = value;

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
	this->formatID[this->formatPointer++] = value;

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

bool BCFEntry::parse(void){
	//std::cerr << this->body->CHROM << ':' << this->body->POS+1 << '\t';
	U32 internal_pos = sizeof(body_type);
	this->__parseID(internal_pos);
	this->__parseRefAlt(internal_pos);
	this->SetRefAlt();

	// start of FORMAT
	this->filter_start = internal_pos;

	// Move up to start of FORMAT
	internal_pos = this->body->l_shared + sizeof(U32)*2;

	// Format key and value
	const base_type& fmt_key = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
#if BCF_ASSERT == 1
	assert(fmt_key.high == 1);
#endif
	this->getInteger(fmt_key.low, internal_pos);
	//this->formatID[this->formatPointer++] = value;

	//std::cerr << (int)fmt_key.high << '\t' << (int)fmt_key.low << '\t' << "id: " << value << std::endl;

	const base_type& fmt_type = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
	//std::cerr << (int)fmt_type.high << '\t' << (int)fmt_type.low << std::endl;

	if(fmt_type.high != 2){
		std::cerr << Helpers::timestamp("LOG","BCF") << "Dropping non-diploid variant: " << this->body->CHROM << ':' << this->body->POS+1 << std::endl;
		this->isGood = false;
		return false;
	}

	// Store virtual offsets into the stream
	// Parse filter

	// Parse info
	// parse format

	this->isGood = true;
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

		if((fmt_type_value1 >> 1) == 0 || (fmt_type_value2 >> 1) == 0) ++n_missing;
	}
	return((double)n_missing/samples);
}

void BCFEntry::SetRefAlt(void){
	this->ref_alt = 0;
	// Set mock ref-alt if not simple
	if(this->alleles[0].length != 1 || this->alleles[1].length != 1){
		this->ref_alt ^= Tomahawk::Constants::REF_ALT_N << 4;
		this->ref_alt ^= Tomahawk::Constants::REF_ALT_N;
		return;
	}

	switch(this->alleles[0].data[0]){
	case 'A': this->ref_alt ^= Tomahawk::Constants::REF_ALT_A << 4; break;
	case 'T': this->ref_alt ^= Tomahawk::Constants::REF_ALT_T << 4; break;
	case 'G': this->ref_alt ^= Tomahawk::Constants::REF_ALT_G << 4; break;
	case 'C': this->ref_alt ^= Tomahawk::Constants::REF_ALT_C << 4; break;
	default:
		std::cerr << Helpers::timestamp("ERROR", "BCF") << "Illegal SNV reference..." << std::endl;
		exit(1);
	}

	switch(this->alleles[1].data[0]){
	case 'A': this->ref_alt ^= Tomahawk::Constants::REF_ALT_A << 0; break;
	case 'T': this->ref_alt ^= Tomahawk::Constants::REF_ALT_T << 0; break;
	case 'G': this->ref_alt ^= Tomahawk::Constants::REF_ALT_G << 0; break;
	case 'C': this->ref_alt ^= Tomahawk::Constants::REF_ALT_C << 0; break;
	case '.': this->ref_alt ^= Tomahawk::Constants::REF_ALT_N << 0; break;
	default:
		std::cerr << Helpers::timestamp("ERROR", "BCF") << "Illegal SNV alt..." << std::endl;
		exit(1);
	}
}

}
}
