#ifndef BCFREADER_H_
#define BCFREADER_H_

#include <cassert>

#define BCF_ASSERT 1

#include "../BasicBuffer.h"
#include "../BGZFController.h"

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

	BCFEntry(void):
		pointer(0),
		limit(262144),
		l_ID(0),
		p_genotypes(0),
		data(new char[this->limit]),
		body(reinterpret_cast<body_type*>(this->data)),
		alleles(new string_type[100]),
		ID(nullptr),
		genotypes(nullptr)
	{

	}

	~BCFEntry(void){ delete [] this->data; }
	void resize(const U32 size){
		char* temp = this->data;
		this->data = new char[size];
		memcpy(this->data, temp, this->pointer);
		std::swap(temp, this->data);
		delete [] temp;
		this->body = reinterpret_cast<body_type*>(this->data);

		if(size > this->limit)
			this->limit = size;
	}

	void add(const char* const data, const U32 length){
		if(this->pointer + length > this-> capacity())
			this->resize(this->pointer + length + 65536);

		memcpy(&this->data[this->pointer], data, length);
		this->pointer += length;
	}

	inline void reset(void){ this->pointer = 0; }
	inline const U32& size(void) const{ return(this->pointer); }
	inline const U32& capacity(void) const{ return(this->limit); }
	U64 sizeBody(void) const{ return(this->body->l_shared + this->body->l_indiv); }

	inline const bool isSimple(void) const{
		return((this->body->n_allele == 2) && (this->alleles[0].length == 1 && this->alleles[1].length == 1));
	}

	void __parseID(U32& internal_pos){
		// Parse ID
		const base_type& ID_base = *reinterpret_cast<const base_type* const>(&this->data[internal_pos]);
		++internal_pos;
		this->ID = &this->data[internal_pos];
		this->l_ID = ID_base.high;
		if(ID_base.high == 0){
			// has no name
			this->l_ID = 0;
		} else if(ID_base.high == 15){
			// next byte is the length array
			std::cerr << "in array" << std::endl;
			const BCFAtomicS32& length = *reinterpret_cast<const BCFAtomicS32* const>(&this->data[internal_pos]);
			internal_pos += sizeof(U32);
			this->l_ID = length.high;
#if BCF_ASSERT == 1
			assert(length.low == 7);
#endif
		}
#if BCF_ASSERT == 1
		assert(ID_base.low == 7);
#endif
		internal_pos += this->l_ID;
	}

	void __parseRefAlt(U32& internal_pos){
		// Parse REF-ALT
		for(U32 i = 0; i < this->body->n_allele; ++i){
			const base_type& alelle_base = *reinterpret_cast<const base_type* const>(&this->data[internal_pos]);
			this->alleles[i].length = alelle_base.high;
			++internal_pos;

			if(alelle_base.low == 15){
				std::cerr << "in array" << std::endl;
				// next byte is the length array
				const BCFAtomicS32& length = *reinterpret_cast<const BCFAtomicS32* const>(&this->data[internal_pos]);
				internal_pos += sizeof(U32);
				this->alleles[i].length = length.high;
#if BCF_ASSERT == 1
				assert(length.low == 7);
#endif
			}
#if BCF_ASSERT == 1
			assert(alelle_base.low == 7);
#endif

			this->alleles[i].data = &this->data[internal_pos];
			internal_pos += this->alleles[i].length;
		}
	}

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

	bool parse(void){
		U32 internal_pos = sizeof(body_type);
		this->__parseID(internal_pos);
		this->__parseRefAlt(internal_pos);
		this->SetRefAlt();

//		if(this->l_ID == 0) std::cerr << '.' << std::endl;
//		else std::cerr << std::string(this->ID, this->l_ID) << std::endl;

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

		// Todo: move out
		this->__parseGenotypes<U32>();
		// Todo: move to RLE parser

		return true;
	}

	void SetRefAlt(void){
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

class BCFReader{
	typedef BCFReader self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::BGZFController bgzf_controller_type;
	typedef IO::BGZFHeader bgzf_type;
	typedef VCF::VCFHeader header_type;
	typedef VCF::VCFHeaderContig contig_type;
	typedef BCFEntry entry_type;

public:
	BCFReader();
	~BCFReader();

	bool basicStats(const entry_type& entry);

	bool nextBlock(void);
	bool nextVariant(BCFEntry& entry);
	bool parseHeader(void);
	bool open(const std::string input);

public:
	std::ifstream stream;
	U64 filesize;
	U32 current_pointer;
	buffer_type buffer;
	buffer_type output_buffer;
	buffer_type header_buffer;
	bgzf_controller_type bgzf_controller;
	header_type header;
};

}
}

#endif /* BCFREADER_H_ */
