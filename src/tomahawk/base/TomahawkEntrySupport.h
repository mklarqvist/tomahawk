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
	typedef struct __support_allele_info{
		U32 n_allele;
		char* allele;
	} allele_info;

public:
	TomahawkSupport(void) :
		virtual_offset(0),
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

public:
	U32 virtual_offset; // virtual offset to end of this value
	float QUAL;
	U16 n_allele;
	// ID
	// byte length of ID
	// even though n_ID is defined as an U32
	// it can be encoded in the stream as a smaller
	// value
	U32 n_ID;
	char* ID;
	// allele info
	allele_info* alleles;
};

}
}

#endif /* TOMAHAWK_BASE_TOMAHAWKENTRYSUPPORT_H_ */
