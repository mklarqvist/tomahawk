#ifndef BCFREADER_H_
#define BCFREADER_H_

#include "../BasicBuffer.h"
#include "../BGZFController.h"

namespace Tomahawk {
namespace IO {

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
		body(reinterpret_cast<const body_type*>(this->data)),
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
		this->body = reinterpret_cast<const body_type*>(this->data);

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
			assert(length.low == 7);
		}
		assert(ID_base.low == 7);
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
				assert(length.low == 7);
			}
			assert(alelle_base.low == 7);

			this->alleles[i].data = &this->data[internal_pos];
			internal_pos += this->alleles[i].length;
		}
	}

	void __parseGenotypes(void){
		U32 internal_pos = this->p_genotypes;
		for(U32 i = 0; i < 32470 * 2; i+=2){
			const SBYTE& fmt_type_value1 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
			const SBYTE& fmt_type_value2 = *reinterpret_cast<SBYTE*>(&this->data[internal_pos++]);
			//std::cerr << (int)((fmt_type_value1>>1) - 1) << ((fmt_type_value2 & 1) == 1 ? '|' : '/') << (int)((fmt_type_value2>>1) - 1) << '\t';
		}
		assert(internal_pos == this->size());
	}

	bool parse(void){
		U32 internal_pos = sizeof(body_type);
		this->__parseID(internal_pos);
		this->__parseRefAlt(internal_pos);

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
		return true;
	}

public:
	U32 pointer; // byte width
	U32 limit;   // capacity
	U32 l_ID;
	U32 p_genotypes; // position genotype data begin

	char* data; // hard copy data to buffer, interpret internally
	const body_type* body; // BCF2 body
	string_type* alleles; // pointer to pointer of ref alleles and their lengths
	char* ID;
	SBYTE* genotypes;
};

class BCFReader{
	typedef BCFReader self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::BGZFController bgzf_controller_type;
	typedef BGZFHeader bgzf_type;
	typedef VCF::VCFHeader header_type;
	typedef VCF::VCFHeaderContig contig_type;

public:
	BCFReader() : filesize(0), current_pointer(0){}
	BCFReader(const std::string file);
	~BCFReader(){}

	bool nextBlock(void){
		// Stream died
		if(!this->stream.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Stream died!" << std::endl;
			return false;
		}

		// EOF
		if(this->stream.tellg() == this->filesize){
			//std::cerr << "eof" << std::endl;
			return false;
		}

		buffer.resize(sizeof(bgzf_type));
		this->stream.read(&buffer.data[0], Constants::BGZF_BLOCK_HEADER_LENGTH);
		const bgzf_type* h = reinterpret_cast<const bgzf_type*>(&buffer.data[0]);
		buffer.pointer = Constants::BGZF_BLOCK_HEADER_LENGTH;
		if(!h->Validate()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Failed to validate!" << std::endl;
			std::cerr << *h << std::endl;
			return false;
		}

		buffer.resize(h->BSIZE); // make sure all data will fit

		// Recast because if buffer is resized then the pointer address is incorrect
		// resulting in segfault
		h = reinterpret_cast<const bgzf_type*>(&buffer.data[0]);

		this->stream.read(&buffer.data[Constants::BGZF_BLOCK_HEADER_LENGTH], (h->BSIZE + 1) - Constants::BGZF_BLOCK_HEADER_LENGTH);
		if(!this->stream.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Truncated file..." << std::endl;
			return false;
		}

		buffer.pointer = h->BSIZE + 1;
		const U32 uncompressed_size = *reinterpret_cast<const U32*>(&buffer[buffer.pointer -  sizeof(U32)]);
		output_buffer.resize(uncompressed_size);
		this->output_buffer.reset();

		if(!this->bgzf_controller.Inflate(buffer, output_buffer)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Failed inflate!" << std::endl;
			return false;
		}

		// BGZF EOF marker
		if(this->output_buffer.size() == 0)
			return false;

		// Reset buffer
		this->buffer.reset();
		this->current_pointer = 0;

		return true;
	}

	bool nextVariant(BCFEntry& entry){
		//BCFEntry entry;
		if(this->current_pointer + 8 > this->output_buffer.size()){
			const U32 partial = this->output_buffer.size() - this->current_pointer;
			entry.add(&this->output_buffer[this->current_pointer], this->output_buffer.size() - this->current_pointer);
			if(!this->nextBlock()){
				std::cerr << "failed to get next block" << std::endl;
				return false;
			}

			entry.add(&this->output_buffer[0], 8 - partial);
			this->current_pointer = 8 - partial;
		} else {
			entry.add(&this->output_buffer[this->current_pointer], 8);
			this->current_pointer += 8;
		}

		U64 remainder = entry.sizeBody();
		//std::cerr << *entry.body << std::endl;
		while(remainder > 0){
			if(this->current_pointer + remainder > this->output_buffer.size()){
				//std::cerr << "partial: " << this->current_pointer + remainder << '/' << this->output_buffer.size() << '\t' << this->output_buffer.size() - this->current_pointer << '\t' << this->current_pointer << std::endl;
				entry.add(&this->output_buffer[this->current_pointer], this->output_buffer.size() - this->current_pointer);
				remainder -= this->output_buffer.size() - this->current_pointer;
				if(!this->nextBlock()){
					std::cerr << "failed to get next block" << std::endl;
					return false;
				}
			} else {
				entry.add(&this->output_buffer[this->current_pointer], remainder);
				this->current_pointer += remainder;
				remainder = 0;
				break;
			}
			//std::cerr << "remaidner now: " << remainder << std::endl;
		}
		//std::cerr << "loaded all: " << remainder << std::endl;
		std::cerr << this->header.getContig(entry.body->CHROM).name << ":" << entry.body->POS << std::endl;
		entry.parse();

		return true;
	}

	bool parseHeader(void){
		if(this->output_buffer.size() == 0){
			std::cerr << "no buffer" << std::endl;
			return false;
		}

		if(strncmp(&this->output_buffer.data[0], "BCF\2\2", 5) != 0){ // weird: should be BCF/2/1
			std::cerr << (int)this->output_buffer[3] << '\t' << (int)this->output_buffer[4] << std::endl;
			std::cerr << "failed to validate" << std::endl;
			return false;
		}

		const U32 l_text = *reinterpret_cast<const U32* const>(&this->output_buffer[5]) + 4;
		this->header_buffer.resize(l_text);

		if(l_text - 5 < this->output_buffer.size()){
			this->header_buffer.Add(&this->output_buffer[5], l_text);
			return true;
		}

		U32 head_read = this->output_buffer.size() - 5;
		this->header_buffer.Add(&this->output_buffer[5], this->output_buffer.size() - 5);

		//U32 p = 0;
		while(this->nextBlock()){
			if(head_read + this->output_buffer.size() >= l_text){
				std::cerr << "remainder: " << l_text - head_read << " and data: " << this->output_buffer.size() << std::endl;
				this->header_buffer.Add(&this->output_buffer[0], l_text - head_read);
				this->current_pointer = l_text - head_read;
				break;
			}
			head_read += this->output_buffer.size();
			this->header_buffer.Add(&this->output_buffer[0], this->output_buffer.size());
		}


		if(!this->header.parse(&this->header_buffer[0], this->header_buffer.size())){
			std::cerr << "failed to parse header" << std::endl;
			return false;
		}

		return true;
	}

	bool open(const std::string input){
		this->stream.open(input, std::ios::binary | std::ios::in | std::ios::ate);
		if(!this->stream.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Failed to open file: " << input << std::endl;
			return false;
		}

		this->filesize = this->stream.tellg();
		this->stream.seekg(0);

		if(!this->stream.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "BCF") << "Bad stream!" << std::endl;
			return false;
		}

		if(!this->nextBlock()){
			std::cerr << "failed ot get first block" << std::endl;
			return false;
		}

		if(!this->parseHeader()){
			std::cerr << "failed to parse bcf header" << std::endl;
			return false;
		}

		return true;
	}

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
