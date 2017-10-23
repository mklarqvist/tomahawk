#ifndef TOMAHAWK_BASE_TOMAHAWKENTRYSUPPORT_H_
#define TOMAHAWK_BASE_TOMAHAWKENTRYSUPPORT_H_

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

	typedef struct __support_allele_info{
		typedef __support_allele_info self_type;
		typedef IO::BasicBuffer buffer_type;

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

		friend buffer_type operator+=(buffer_type& buffer, const self_type entry){
			BYTE out = entry.length << 2;
			out |= entry.type;
			buffer += out;
			return(buffer);
		}

		BYTE type: 2, length: 6;
	} typed_value;

public:
	explicit TomahawkSupport(void) :
		l_body(0),
		QUAL(0),
		n_allele(0),
		n_ID(0),
		ID(nullptr),
		alleles(nullptr)
	{}

	~TomahawkSupport(void){
		delete [] this->alleles;
	}

	bool parse(void);
	bool parseID(void);
	bool parseAlleles(void);

	// Write out entry using BCF entry as template
	// and injects into buffer
	bool write(const bcf_type& entry, buffer_type& buffer){
		// Determine offset
		// Base length
		// offset + QUAL + n_alleles
		U32 offset = sizeof(U32) + sizeof(float) + sizeof(U16);
		const U32 buffer_start = buffer.pointer;

		// ID length
		if(entry.l_ID < 63) offset += sizeof(BYTE);
		else if(entry.l_ID < 256) offset += 2*sizeof(BYTE); // BYTE + BYTE
		else offset += sizeof(BYTE) + sizeof(U16); // BYTE + U16

		// Allele length
		for(U32 i = 0; i < entry.body->n_allele; ++i){
			if(entry.alleles[i].length < 63) offset += sizeof(BYTE);
			else if(entry.alleles[i].length < 256) offset += 2*sizeof(BYTE); // BYTE + BYTE
			else offset += sizeof(BYTE) + sizeof(U16); // BYTE + U16
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
		//std::cerr << (S32)entry.l_ID << std::endl;
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

		// Write out alleles
		for(U32 i = 0; i < entry.body->n_allele; ++i){
			//std::cerr << (S32)entry.alleles[i].length << std::endl;
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
		}

		//std::cerr << "Expected: " << offset << "; observed: " << (S32)buffer.pointer - offset << std::endl;
		assert((S32)buffer.pointer - (buffer_start + offset) == 0);

		return true;
	}

public:
	// relative virtual offset to end of this value
	// this is equivalent to its length
	U32 l_body;
	float QUAL;
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
