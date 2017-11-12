#ifndef BCF_BCFENTRY_H_
#define BCF_BCFENTRY_H_

#include "../vcf/VCFHeader.h"
#include "../../third_party/xxhash/xxhash.h"
#include "../BasicBuffer.h"

namespace Tomahawk {
namespace BCF {

#define BCF_ASSERT 1
#define BCF_HASH_SEED 452930477

const BYTE BCF_UNPACK_TOMAHAWK[3] = {2, 0, 1};
#define BCF_UNPACK_GENOTYPE(A) BCF_UNPACK_TOMAHAWK[(A >> 1)]
const char BCF_TYPE_SIZE[8] = {0,1,2,4,0,4,0,1};

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
	U32 n_sample: 24, n_fmt: 8;
};

struct BCFTypeString{
	U16 length;
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

	inline void reset(void){
		this->pointer = 0;
		this->isGood = false;
		this->infoPointer = 0;
		this->formatPointer = 0;
		this->filterPointer = 0;
		this->n_filter = 0;
		this->filter_start = 0;
	}

	inline const U32& size(void) const{ return(this->pointer); }
	inline const U32& capacity(void) const{ return(this->limit); }
	inline const U64 sizeBody(void) const{ return(this->body->l_shared + this->body->l_indiv); }

	inline const bool isBiallelicSimple(void) const{
		return((this->body->n_allele == 2) && (this->alleles[0].length == 1 && this->alleles[1].length == 1));
	}
	inline const bool isBiallelic(void) const{ return(this->body->n_allele == 2); }
	inline const bool isSimple(void) const{
		return(this->alleles[0].length == 1 && this->alleles[1].length == 1);
	}

	void __parseID(U32& internal_pos);
	void __parseRefAlt(U32& internal_pos);

	bool parse(void);
	void SetRefAlt(void);
	double getMissingness(const U64& samples) const;
	inline const bool& good(void) const{ return(this->isGood); }

	inline const SBYTE& getSBYTE(U32& pos){ return(*reinterpret_cast<const SBYTE* const>(&this->data[pos++])); }

	inline const S16& getS16(U32& pos){
		const S16& t = *reinterpret_cast<const S16* const>(&this->data[pos]);
		pos += sizeof(S16);
		return(t);
	}

	inline const S32& getS32(U32& pos){
		const S32& t = *reinterpret_cast<const S32* const>(&this->data[pos]);
		pos += sizeof(S32);
		return(t);
	}

	inline const S32 getInteger(const BYTE& key, U32& pos){
		S32 value = 0;
		switch(key){
		case(1): value = *reinterpret_cast<const SBYTE* const>(&this->data[pos++]); break;
		case(2): value = *reinterpret_cast<const S16* const>(&this->data[pos]); pos+=sizeof(S16);  break;
		case(3): value = *reinterpret_cast<const S32* const>(&this->data[pos]); pos+=sizeof(S32);  break;
		case(0): std::cerr << "FLAG value" << std::endl; break;
		default: std::cerr << "illegal" << std::endl; exit(1);
		}
		return value;
	}

	inline const float getFloat(U32& pos){
		const float val = *reinterpret_cast<const float* const>(&this->data[pos]);
		pos += sizeof(float);
		return val;
	}

	inline const char getChar(U32& pos){
		return(*reinterpret_cast<const char* const>(&this->data[pos++]));
	}

	inline const U64 hashFilter(void){return(XXH64((const void*)this->filterID, sizeof(U32)*this->filterPointer, BCF_HASH_SEED));}
	inline const U64 hashInfo(void){return(XXH64((const void*)this->infoID, sizeof(U32)*this->infoPointer, BCF_HASH_SEED));}
	inline const U64 hashFormat(void){return(XXH64((const void*)this->formatID, sizeof(U32)*this->formatPointer, BCF_HASH_SEED));}

	// Iterators over fields
	bool nextFilter(S32& value, U32& position);
	bool nextInfo(S32& value, U32& length, BYTE& value_type, U32& position);
	bool nextFormat(S32& value, U32& length, BYTE& value_type, U32& position);

public:
	U32 pointer; // byte width
	U32 limit;   // capacity
	U32 l_ID;
	U32 p_genotypes; // position genotype data begin
	BYTE ref_alt; // parsed
	bool isGood;
	char* data; // hard copy data to buffer, interpret internally
	body_type* body; // BCF2 body
	string_type* alleles; // pointer to pointer of ref alleles and their lengths
	char* ID;
	SBYTE* genotypes;

	//
	U32 filter_start;
	U32 n_filter;
	base_type filter_key;

	// Vectors of identifiers
	U16 filterPointer;
	U16 infoPointer;
	U16 formatPointer;
	// FILTER
	U32* filterID;
	// INFO
	U32* infoID;
	// FORMAT
	U32* formatID;
};

}
}

#endif /* BCF_BCFENTRY_H_ */
