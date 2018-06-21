#ifndef BCFENTRY_H_
#define BCFENTRY_H_

#include "io/basic_buffer.h"

namespace tomahawk {
namespace bcf {

const BYTE BCF_UNPACK_TOMAHAWK[3] = {2, 0, 1};
#define BCF_UNPACK_GENOTYPE(A) BCF_UNPACK_TOMAHAWK[(A >> 1)]
const char BCF_TYPE_SIZE[8] = {0,1,2,4,0,4,0,1};

#pragma pack(push, 1)
struct __attribute__((packed, aligned(1))) BCFAtomicBase{
	BYTE low: 4, high: 4;
};

struct __attribute__((packed, aligned(1))) BCFEntryBody{
	typedef BCFEntryBody self_type;

	BCFEntryBody(); // disallow ctor
	~BCFEntryBody(); // disallow dtor

	// For debugging only
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
#pragma pack(pop)

struct BCFTypeString{
	typedef BCFTypeString self_type;

	BCFTypeString() : length(0), data(nullptr){}
	BCFTypeString(const self_type& other) : length(other.length), data(new char[other.length]){ memcpy(this->data, other.data, other.length); }
	void operator=(const self_type& other){
		this->length = other.length;
		delete [] this->data;
		this->data = new char[other.length];
		memcpy(this->data, other.data, other.length);
	}
	void operator()(const char* const data, const U16 length){
		this->length = length;
		delete [] this->data;
		this->data = new char[length];
		memcpy(this->data, data, length);
	}
	~BCFTypeString(){
		delete [] this->data;
	}

	U16 length;
	char* data;
};

struct BCFKeyTuple{
	BCFKeyTuple() : mapID(-1), primitive_type(0), l_stride(0), l_offset(0){}
	~BCFKeyTuple(){}

	S32  mapID;
	BYTE primitive_type;
	U32  l_stride;
	U32  l_offset;
};

struct BCFGenotypeSupport{
private:
	typedef BCFGenotypeSupport self_type;

public:
	BCFGenotypeSupport() :
		hasGenotypes(false),
		hasMissing(false),
		hasEOV(false),
		mixedPhasing(false),
		ploidy(0),
		phase(0),
		n_missing(0),
		n_eov(0)
	{}

	BCFGenotypeSupport(const self_type& other) :
		hasGenotypes(other.hasGenotypes),
		hasMissing(other.hasMissing),
		hasEOV(other.hasEOV),
		mixedPhasing(other.mixedPhasing),
		ploidy(other.ploidy),
		phase(other.phase),
		n_missing(other.n_missing),
		n_eov(other.n_eov)
	{
	}

	void operator=(const self_type& other){
		this->hasGenotypes = other.hasGenotypes;
		this->hasMissing   = other.hasMissing;
		this->hasEOV       = other.hasEOV;
		this->mixedPhasing = other.mixedPhasing;
		this->ploidy       = other.ploidy;
		this->phase        = other.phase;
		this->n_missing    = other.n_missing;
		this->n_eov        = other.n_eov;
	}

	void reset(void){
		this->hasGenotypes = false;
		this->hasMissing = false;
		this->hasEOV = false;
		this->mixedPhasing = false;
		this->ploidy = 0;
		this->phase = 0;
		this->n_missing = 0;
		this->n_eov = 0;
	}

public:
	bool hasGenotypes;
	bool hasMissing;
	bool hasEOV;
	bool mixedPhasing;
	U32  ploidy;
	BYTE phase;
	U32  n_missing;
	U32  n_eov;
};

struct BCFEntry{
private:
	typedef BCFEntry           self_type;
	typedef io::BasicBuffer    buffer_type;
	typedef BCFEntryBody       body_type;
	typedef BCFTypeString      string_type;
	typedef BCFAtomicBase      base_type;
	typedef BCFGenotypeSupport gt_support_type;

public:
	BCFEntry(void);  // ctor
	BCFEntry(const self_type& other);  // copy ctor
	BCFEntry(self_type&& other) noexcept;
	BCFEntry& operator=(const self_type& other);
	BCFEntry& operator=(self_type&& other) noexcept;
	~BCFEntry(void); // dtor

	void resize(const U32 size);
	void add(const char* const data, const U32 length);

	inline void reset(void){
		this->l_data        = 0;
		this->isGood        = false;
		this->infoPointer   = 0;
		this->formatPointer = 0;
		this->filterPointer = 0;
		this->n_filter      = 0;
		this->filter_start  = 0;
		this->gt_support.reset();
	}

	inline const U32& size(void) const{ return(this->l_data); }
	inline const U32& capacity(void) const{ return(this->l_capacity); }
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

	bool parse(const U64 n_samples);

	void SetRefAlt(void);

	//
	bool assessGenotypes(const U64 n_samples){
		//if(this->hasGenotypes == false)
		//	return false;

		const BYTE ploidy = this->formatID[0].l_stride;
		// Todo: other primitives
		U64 current_sample = 0;
		const char* const internal_data = &this->data[this->formatID[0].l_offset];
		U32 internal_data_offset = 0;
		this->gt_support.ploidy = ploidy;
		this->gt_support.hasGenotypes = true;

		BYTE first_phase = 0;

		if(this->formatID[0].primitive_type == 1){
			for(U32 p = 0; p < ploidy; ++p){
				const BYTE& ref  = *reinterpret_cast<const BYTE* const>(&internal_data[internal_data_offset]);
				if((ref >> 1) == 0){
					this->gt_support.hasMissing = true;
					++this->gt_support.n_missing;
				} else if(ref == 0x81){
					this->gt_support.hasEOV = true;
					++this->gt_support.n_eov;
				}
				if(p + 1 == ploidy) first_phase = ref & 1;
				internal_data_offset += sizeof(BYTE);
			}
			++current_sample;

			for(U32 i = ploidy; i < n_samples*ploidy; i+=ploidy, ++current_sample){
				// retrieve ploidy primitives
				for(U32 p = 0; p < ploidy; ++p){
					const BYTE& ref  = *reinterpret_cast<const BYTE* const>(&internal_data[internal_data_offset]);
					if((ref >> 1) == 0){
						//std::cerr << "is missing" << std::endl;
						this->gt_support.hasMissing = true;
						++this->gt_support.n_missing;
					} else if(ref == 0x81){
						//std::cerr << "is vector eof" << std::endl;
						this->gt_support.hasEOV = true;
						++this->gt_support.n_eov;
					}

					if(p + 1 == ploidy){
						if(first_phase != (ref & 1)){
							//std::cerr << "triggering mixed phase" << std::endl;
							this->gt_support.mixedPhasing = true;
							this->gt_support.phase = 0;
						}
					}
					internal_data_offset += sizeof(BYTE);
				}
			}

			if(this->gt_support.mixedPhasing == false){
				this->gt_support.phase = first_phase;
			}
		}
		return true;
	}

	inline const bool& good(void) const{ return(this->isGood); }

	/**<
	 * Decode an integer primitive from a BCF buffer stream. Forces all
	 * return types to be of type S32 and with missing and EOV values
	 * expanded to match this possibly larger primitive type.
	 * @param key
	 * @param pos
	 * @return
	 */
	inline const S32 getInteger(const BYTE& key, U32& pos){
		S32 value = 0;
		if(key == 1){
			const SBYTE& ref  = *reinterpret_cast<const SBYTE* const>(&this->data[pos++]);
			const BYTE&  uref = *reinterpret_cast<const BYTE* const>(&ref);
			if(uref == 0x80){
				//std::cerr << "is missing" << std::endl;
				return(value = 0x80000000);
			} else if(uref == 0x81){
				return(value = 0x80000001);
				//std::cerr << "is vector eof" << std::endl;
			}
			return(value = ref);
		} else if(key == 2){
			const S16& ref  = *reinterpret_cast<const S16* const>(&this->data[pos]);
			const U16& uref = *reinterpret_cast<const U16* const>(&ref);
			pos+=sizeof(S16);

			if(uref == 0x8000){
				//std::cerr << "is missing s16" << std::endl;
				return(value = 0x80000000);
			} else if(uref == 0x8001){
				//std::cerr << "is vector eof" << std::endl;
				return(value = 0x80000001);
			}
			return(value = ref);
		} else if(key == 3){
			const S32& ref  = *reinterpret_cast<const S32* const>(&this->data[pos]);
			const U32& uref = *reinterpret_cast<const U32* const>(&ref);
			pos+=sizeof(S32);

			if(uref == 0x80000000){
				//std::cerr << "is missing" << std::endl;
				return(value = 0x80000000);
			} else if(uref == 0x80000001){
				//std::cerr << "is vector eof" << std::endl;
				return(value = 0x80000001);
			}
			return(value = ref);
		} else if(key == 0){
			return 0;
		} else {
			std::cerr << "illegal type" << std::endl;
			exit(1);
		}
	}

	/**<
	 * Decodes a float values from a BCF buffer stream
	 * @param pos
	 * @return
	 */
	inline const float getFloat(U32& pos){
		const float& val = *reinterpret_cast<const float* const>(&this->data[pos]);
		pos += sizeof(float);
		return val;
	}

	/**<
	 * Decodes a char from a BCF buffer stream
	 * @param pos
	 * @return
	 */
	inline const char getChar(U32& pos){ return(this->data[pos++]); }
	inline const char* const getCharPointer(U32& pos){ return(&this->data[pos]); }

	// Iterators over fields
	bool nextFilter(S32& value, U32& position);
	bool nextInfo(S32& value, U32& length, BYTE& value_type, U32& position);
	bool nextFormat(S32& value, U32& length, BYTE& value_type, U32& position);


public:
	U32 l_data;     // byte width
	U32 l_capacity;       // capacity
	U32 l_ID;
	BYTE ref_alt;    // parsed
	bool isGood;
	char* data;      // hard copy data to buffer, interpret internally
	body_type* body; // BCF2 body
	string_type* alleles; // pointer to pointer of ref alleles and their lengths
	char* ID;

	bool hasGenotypes;
	BYTE ploidy;

	//
	U32 filter_start;
	U32 n_filter;
	base_type filter_key;

	// Vectors of identifiers
	U16 filterPointer;
	U16 infoPointer;
	U16 formatPointer;
	gt_support_type gt_support;

	// FILTER
	BCFKeyTuple* filterID;
	BCFKeyTuple* infoID;
	BCFKeyTuple* formatID;
};

}
}

#endif /* BCFENTRY_H_ */
