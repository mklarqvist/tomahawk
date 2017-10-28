#include "TomahawkEntrySupport.h"

namespace Tomahawk{
namespace Support{

TomahawkSupport::TomahawkSupport(void) :
		l_body(0),
		QUAL(0),
		n_allele(0),
		n_ID(0),
		ID(nullptr),
		alleles(nullptr)
	{}

TomahawkSupport::~TomahawkSupport(void){
	delete [] this->alleles;
}

bool TomahawkSupport::write(const bcf_type& entry, buffer_type& buffer){
	// Determine offset
	// Base length
	// offset + QUAL + n_alleles
	U32 offset = sizeof(U32) + sizeof(float) + sizeof(U16);
	const U32 buffer_start = buffer.pointer;

	// ID length
	if(entry.l_ID < 63) offset += sizeof(BYTE);
	else if(entry.l_ID < 256) offset += 2*sizeof(BYTE); // BYTE + BYTE
	else offset += sizeof(BYTE) + sizeof(U16); // BYTE + U16
	offset += entry.l_ID;

	// Allele length
	for(U32 i = 0; i < entry.body->n_allele; ++i){
		if(entry.alleles[i].length < 63) offset += sizeof(BYTE);
		else if(entry.alleles[i].length < 256) offset += 2*sizeof(BYTE); // BYTE + BYTE
		else offset += sizeof(BYTE) + sizeof(U16); // BYTE + U16
		offset += entry.alleles[i].length;
	}

	// Assert that data will fit in buffer
	if(buffer.pointer + offset > buffer.width)
		buffer.resize(buffer.width * 1.2);

	// Write out data
	// offset is
	buffer += offset;
	buffer += entry.body->QUAL;
	buffer += (U16)entry.body->n_allele;

	// Write out ID
	typed_value n_ID;
	if(entry.l_ID < 63){
		n_ID.type = typed_value::BYTE_TYPE;
		n_ID.length = entry.l_ID;
		buffer += n_ID;
	} else if(entry.l_ID < 256){
		n_ID.type = typed_value::BYTE_TYPE;
		n_ID.length = 63;
		buffer += n_ID;
		buffer += (BYTE)entry.l_ID;
	} else{
		n_ID.type = typed_value::U16_TYPE;
		n_ID.length = 63;
		buffer += n_ID;
		buffer += (U16)entry.l_ID;
	}
	buffer.Add(entry.ID, entry.l_ID);

	// Write out alleles
	for(U32 i = 0; i < entry.body->n_allele; ++i){
		// Write out allele
		typed_value n_ID;
		if(entry.alleles[i].length < 63){
			n_ID.type = typed_value::BYTE_TYPE;
			n_ID.length = entry.alleles[i].length;
			buffer += n_ID;
		} else if(entry.alleles[i].length < 256){
			n_ID.type = typed_value::BYTE_TYPE;
			n_ID.length = 63;
			buffer += n_ID;
			buffer += (BYTE)entry.alleles[i].length;
		} else{
			n_ID.type = typed_value::U16_TYPE;
			n_ID.length = 63;
			buffer += n_ID;
			buffer += (U16)entry.alleles[i].length;
		}
		buffer.Add(entry.alleles[i].data, entry.alleles[i].length);
	}

	assert((S32)buffer.pointer - (buffer_start + offset) == 0);

	return true;
}

}
}
