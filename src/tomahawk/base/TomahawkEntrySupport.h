#ifndef TOMAHAWK_BASE_TOMAHAWKENTRYSUPPORT_H_
#define TOMAHAWK_BASE_TOMAHAWKENTRYSUPPORT_H_

#include <cassert>

#include "../../io/BasicBuffer.h"
#include "../../io/bcf/BCFEntry.h"

namespace Tomahawk{
namespace Support{

// Do NOT reinterpret_cast this struct as an array
// as offsets needs to be interpreted
#pragma pack(1)
struct TomahawkSupport{
private:
	typedef TomahawkSupport self_type;
	typedef BCF::BCFEntry bcf_type;
	typedef IO::BasicBuffer buffer_type;

public:
	typedef struct __support_allele_info{
		typedef __support_allele_info self_type;
		typedef IO::BasicBuffer buffer_type;

		explicit __support_allele_info(void) : l_allele(0), allele(nullptr){}
		~__support_allele_info(){
			// do not clear char pointer
			// points to byte stream
		}

		friend buffer_type& operator+=(buffer_type& buffer, const self_type& entry){
			// Write out allele
			typed_value n_ID;
			if(entry.l_allele < 63){
				n_ID.type = typed_value::BYTE_TYPE;
				n_ID.length = entry.l_allele;
				buffer += n_ID;
			}
			else if(entry.l_allele < 256){
				n_ID.type = typed_value::BYTE_TYPE;
				n_ID.length = 63;
				buffer += n_ID;
				buffer += (BYTE)entry.l_allele;
			}
			else{
				n_ID.type = typed_value::U16_TYPE;
				n_ID.length = 63;
				buffer += n_ID;
				buffer += (U16)entry.l_allele;
			}

			return(buffer);
		}

		U16 l_allele;
		char* allele;
	} allele_info;

	typedef struct __typedValue{
		typedef __typedValue self_type;
		typedef IO::BasicBuffer buffer_type;
		enum typedValueType{BYTE_TYPE, U16_TYPE, U32_TYPE, U64_TYPE};

		explicit __typedValue(void): type(0), length(0){}
		~__typedValue(){}

		inline const bool isEmpty(void) const{ return(true); }
		inline const bool isFull(void) const{ return(this->length == 63); }

		friend buffer_type& operator+=(buffer_type& buffer, const self_type& entry){
			BYTE out = entry.length << 2;
			out |= entry.type;
			buffer += out;
			return(buffer);
		}

		BYTE type: 2, length: 6;
	} typed_value;

public:
	explicit TomahawkSupport(void);
	~TomahawkSupport(void);

	// Parse everything in this entry
	bool parse(void);

	// Parse the name only
	bool parseID(void);

	// Parse the alleles only
	bool parseAlleles(void);

	// Write out entry using BCF entry as template
	// and injects into buffer
	bool write(const bcf_type& entry, buffer_type& buffer);

public:
	// Relative virtual offset to end of this
	// (parsed) struct in the byte stream
	// This is equivalent to its length
	U32 l_body;
	float QUAL;
	// Number of alleles
	U16 n_allele;
	// ID
	// byte length of ID
	// even though n_ID is defined as an U16
	// it can be encoded in the stream as a smaller
	// value
	// Names are limited to 16 bits
	U16 n_ID;
	char* ID;
	// allele info
	// ALTs are limited to 16 bits each
	allele_info* alleles;
};

}
}

#endif /* TOMAHAWK_BASE_TOMAHAWKENTRYSUPPORT_H_ */
