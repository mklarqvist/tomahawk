#ifndef BCFENTRY_H_
#define BCFENTRY_H_

#include "../BasicBuffer.h"

namespace Tomahawk {
namespace BCF {

const BYTE BCF_UNPACK_TOMAHAWK[3] = {2, 0, 1};
#define BCF_UNPACK_GENOTYPE(A) BCF_UNPACK_TOMAHAWK[(A >> 1)]


#pragma pack(1)
struct BCFAtomicBase{
	BYTE low: 4, high: 4;
};

#pragma pack(1)
struct BCFAtomicSBYTE{
	SBYTE low: 4, high: 4;
};

#pragma pack(1)
struct BCFAtomicS16{
	S16 low: 4, high: 12;
};

#pragma pack(1)
struct BCFAtomicS32{
	S32 low: 4, high: 28;
};

#pragma pack(1)
struct BCFEntryBody{
	typedef BCFEntryBody self_type;

	BCFEntryBody(); // disallow ctor and dtor
	~BCFEntryBody();

	friend std::ostream& operator<<(std::ostream& os, const self_type& header){
		os << "l_shared\t" << (U32)header.l_shared << '\n';
		os << "l_indiv\t" << (U32)header.l_indiv << '\n';
		os << "CHROM\t" << (U32)header.CHROM << '\n';
		os << "POS\t" << (U32)header.POS << '\n';
		os << "rlen\t" << (S32)header.rlen << '\n';
		os << "QUAL\t" << (U32)header.QUAL << '\n';
		os << "n_allele\t" << (U32)header.n_allele << '\n';
		os << "n_info\t" << (U16)header.n_info << '\n';
		os << "n_fmt\t" << (U32)header.n_fmt << '\n';
		os << "n_sample\t" << (U32)header.n_sample;

		return os;
	}

	U32 l_shared;
	U32 l_indiv;
	S32 CHROM;
	S32 POS;
	S32 rlen;
	float QUAL;
	U32 n_info: 16, n_allele: 16;
	U32 n_sample: 8, n_fmt: 24;
};

struct BCFTypeString{
	typedef BCFAtomicBase base_type;

	U32 length;
	char* data; // reinterpret me as char*
};

struct BCFEntry{
	typedef IO::BasicBuffer buffer_type;
	typedef BCFEntryBody body_type;
	typedef BCFTypeString string_type;
	typedef BCFAtomicBase base_type;

	BCFEntry(void);
	~BCFEntry(void);

	void resize(const U32 size);
	void add(const char* const data, const U32 length);
	inline void reset(void){ this->pointer = 0; }
	inline const U32& size(void) const{ return(this->pointer); }
	inline const U32& capacity(void) const{ return(this->limit); }
	inline U64 sizeBody(void) const{ return(this->body->l_shared + this->body->l_indiv); }

	inline const bool isSimple(void) const{
		return((this->body->n_allele == 2) && (this->alleles[0].length == 1 && this->alleles[1].length == 1));
	}

	void __parseID(U32& internal_pos);
	void __parseRefAlt(U32& internal_pos);

	template <class T>
	void __parseGenotypes(void){
		U32 internal_pos = this->p_genotypes;
		T length = 1;
#if BCF_ASSERT == 1
		U32 sumLength = 0;
#endif
		U32 runs = 0;
		const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
		BYTE packed = (BCF_UNPACK_GENOTYPE(fmt_type_value1) << 2) | BCF_UNPACK_GENOTYPE(fmt_type_value2);

		for(U32 i = 2; i < 32470 * 2; i+=2){
			const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
			const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
			const BYTE packed_internal = (BCF_UNPACK_GENOTYPE(fmt_type_value1) << 2) | BCF_UNPACK_GENOTYPE(fmt_type_value2);
			if(packed != packed_internal){
#if BCF_ASSERT == 1
				sumLength += length;
#endif
				length = 1;
				packed = packed_internal;
				++runs;
				continue;
			}
			++length;
		}
		++runs;
#if BCF_ASSERT == 1
		sumLength += length;
		assert(sumLength == 32470);
		assert(internal_pos == this->size());
#endif
	}

	bool parse(void);
	void SetRefAlt(void);

public:
	U32 pointer; // byte width
	U32 limit;   // capacity
	U32 l_ID;
	U32 p_genotypes; // position genotype data begin
	BYTE ref_alt; // parsed

	char* data; // hard copy data to buffer, interpret internally
	body_type* body; // BCF2 body
	string_type* alleles; // pointer to pointer of ref alleles and their lengths
	char* ID;
	SBYTE* genotypes;
};

}
}

#endif /* BCFENTRY_H_ */
