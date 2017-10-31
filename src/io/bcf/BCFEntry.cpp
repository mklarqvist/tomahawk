#ifndef BCFENTRY_CPP_
#define BCFENTRY_CPP_

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
		this->l_ID = this->getInteger(array_base, this->data, internal_pos);
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
			this->alleles[i].length = this->getInteger(array_base, this->data, internal_pos);
			this->alleles[i].data = &this->data[internal_pos];
		}

		//std::cerr << std::string(this->alleles[i].data, this->alleles[i].length) << '\t';
		internal_pos += this->alleles[i].length;
	}
}

bool BCFEntry::parse(void){
	//std::cerr << this->body->CHROM << ':' << this->body->POS+1 << '\t';
	U32 internal_pos = sizeof(body_type);
	this->__parseID(internal_pos);
	this->__parseRefAlt(internal_pos);
	this->SetRefAlt();

	//std::cerr << this->body->CHROM << ':' << this->body->POS+1 << '\t' << this->body->n_allele << '\t' << this->body->n_fmt << '\t' << this->body->n_info << '\t' << this->body->n_sample << std::endl;

	// At FILTER
	const base_type& filter_key = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
	U32 n_filter = filter_key.high;
	if(n_filter == 15) n_filter = this->getInteger(filter_key, this->data, internal_pos);

	for(U32 i = 0; i < n_filter; ++i){
		const S32 filter_value = this->getInteger(filter_key, this->data, internal_pos);
		this->filterID[this->filterPointer++] = filter_value;
		//std::cerr << "FILTER: " << (int)n_filter << '\t' << (int)filter_key.low << '\t' << "id: " << filter_value << std::endl;
	}

	// At INFO

	for(U32 i = 0; i < this->body->n_info; ++i){
		const base_type& info_key = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
#if BCF_ASSERT == 1
		// This first bit returns a single identifier
		// to a field. It should always be a single
		// value
		assert(info_key.high == 1);
#endif
		const S32 id = this->getInteger(info_key, this->data, internal_pos);
		//std::cerr << i << " INFO: " << (int)info_key.high << '/' << (int)info_key.low << " id: " << id << std::endl;
		//std::cerr << id << ":";

		this->infoID[this->infoPointer++] = id;

		const base_type& info_value = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
		S32 n_info_value = info_value.high;
		if(n_info_value == 15){
			//std::cerr << "is array" << std::endl;
			const base_type& array_base = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
			n_info_value = this->getInteger(array_base, this->data, internal_pos);
			//std::cerr << "width now: " << n_info_value << std::endl;
		}

		//std::cerr << '\t' << (int)info_value.low << '\t' << n_info_value << std::endl;
		if(info_value.low <= 3){
			for(U32 j = 0; j < n_info_value; ++j){
				const S32 val = this->getInteger(info_value, this->data, internal_pos);
				//std::cerr << '\t' << j << ": " << val << std::endl;
				//std::cerr << val << ';';
			}

		} else if(info_value.low == 5){
			for(U32 j = 0; j < n_info_value; ++j){
				//std::cerr << '\t' << j << ": " << this->getFloat(this->data, internal_pos) << std::endl;
				//std::cerr << this->getFloat(this->data, internal_pos) << ';';
				this->getFloat(this->data, internal_pos);
			}

		} else if(info_value.low == 7){
			//std::cerr << '\t';
			for(U32 j = 0; j < n_info_value; ++j){
				//std::cerr << this->getChar(this->data, internal_pos);
				this->getChar(this->data, internal_pos);
			}
			//std::cerr << ';';

		} else {
			std::cerr << "impossible: " << (int)info_value.low << std::endl;
			exit(1);
		}

		//internal_pos += n_info * BCF_TYPE_SIZE[info_key.low];
		//std::cerr << i << " INFO: " << (int)n_info << '\t' << (int)info_key.low << " id: " << 1 << "\tjump: " << 1 << std::endl;
	}
	//std::cerr << std::endl;

	//std::cerr << this->filterPointer << '\t' << this->infoPointer << '\t' << this->formatPointer << std::endl;
	std::cerr << this->hashFilter() << '\t' << this->hashInfo() << '\t' << this->hashFormat() << std::endl;

#if BCF_ASSERT == 1
	// Assert all FILTER and INFO data have been successfully
	// parsed. This is true when the byte pointer equals the
	// start position of the FORMAT fields which are encoded
	// in the meta header structure
	assert(internal_pos == (this->body->l_shared + sizeof(U32)*2));
#endif

	//if(this->body->POS+1 == 118498) exit(1); */

	// Move up to start of FORMAT
	internal_pos = this->body->l_shared + sizeof(U32)*2;

	// Format key and value
	const base_type& fmt_key = *reinterpret_cast<const base_type* const>(&this->data[internal_pos++]);
#if BCF_ASSERT == 1
	assert(fmt_key.high == 1);
#endif
	const S32 value = this->getInteger(fmt_key, this->data, internal_pos);


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
		//std::cerr << i << ':' << " " << (int)fmt_type_value1 << ',' << (int)fmt_type_value2 << '\t' << (int)(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1)) << ',' << (int)(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2)) << std::endl;

		if(fmt_type_value1 < 0 || fmt_type_value2 < 0)
			return(1);

		if(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) == 2 || BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) == 2)
		std::cerr << (int)BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) << '\t' << (int)BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) << std::endl;

		if(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) == 2 || BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) == 2) ++n_missing;
	}
	return((double)n_missing/samples);
}

void BCFEntry::SetRefAlt(void){
	this->ref_alt = 0;
	if(this->alleles[0].length != 1 || this->alleles[1].length != 1){
		//std::cerr << "setting mock refalt for: " << std::string(this->alleles[0].data, this->alleles[0].length) << '\t' << std::string(this->alleles[1].data, this->alleles[1].length) << std::endl;
		this->ref_alt ^= Tomahawk::Constants::REF_ALT_N << 4;
		this->ref_alt ^= Tomahawk::Constants::REF_ALT_N;
		return;
	}

	switch(this->alleles[0].data[0]){
	case 'A': this->ref_alt ^= Tomahawk::Constants::REF_ALT_A << 4; break;
	case 'T': this->ref_alt ^= Tomahawk::Constants::REF_ALT_T << 4; break;
	case 'G': this->ref_alt ^= Tomahawk::Constants::REF_ALT_G << 4; break;
	case 'C': this->ref_alt ^= Tomahawk::Constants::REF_ALT_C << 4; break;
	case '.': this->ref_alt ^= Tomahawk::Constants::REF_ALT_N << 4; break;
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



#endif /* BCFENTRY_CPP_ */
