#ifndef TOMAHAWK_TomahawkBlockManager_H_
#define TOMAHAWK_TomahawkBlockManager_H_

#include "../support/simd_definitions.h"
#include "../algorithm/GenotypeBitPacker.h"
#include "base/TomahawkSupport.h"
#include "../tomahawk/base/TomahawkEntryMeta.h"
#include "../totempole/TotempoleEntry.h"
#include "../totempole/TotempoleReader.h"
#include "base/genotype_container_bitvector.h"

namespace Tomahawk{

template <class T, class Y>
struct TomahawkBlock{
	typedef Y		  	  type;
	typedef type          value_type;
	typedef type         *pointer;
	typedef const type   *const_pointer;
	typedef type         &reference;
	typedef const type   &const_reference;
	typedef size_t        size_type;
	typedef ptrdiff_t     difference_type;

	typedef TomahawkEntryMeta<T> meta_type;

public:
	TomahawkBlock(const char* target, const Totempole::TotempoleEntry& support) :
		metaPointer(0),
		runsPointer(0),
		support(&support),
		meta(reinterpret_cast<const TomahawkEntryMeta<T>* const>(target)),
		runs(reinterpret_cast<const type* const>(&target[(TOMAHAWK_ENTRY_META_SIZE + sizeof(T)) * support.variants])),
		packed(new Base::GenotypeContainerBitvector)
	{

	}

	~TomahawkBlock() noexcept{}

	// copy constructor
	TomahawkBlock(const TomahawkBlock& other) :
		metaPointer(other.metaPointer),
		runsPointer(other.runsPointer),
		support(other.support),
		meta(other.meta),
		runs(other.runs),
		packed(other.packed)
	{

	}

	void operator=(const TomahawkBlock& other){
		this->metaPointer = other.metaPointer;
		this->runsPointer = other.runsPointer;
	}

	// move constructor
	TomahawkBlock(TomahawkBlock&& other) noexcept :
		metaPointer(other.metaPointer),
		runsPointer(other.runsPointer),
		support(other.support),
		meta(other.meta),
		runs(other.runs),
		packed(other.packed)
	{

	}

	inline void updatePacked(const TomahawkBlock& self){
		this->packed = new Base::GenotypeContainerBitvector(self.packed);
	}


	inline void operator++(void){
		this->runsPointer += this->meta[this->metaPointer].runs;
		++this->metaPointer;
	}

	inline void operator--(void){
		--this->metaPointer;
		this->runsPointer -= this->meta[this->metaPointer].runs;
	}

	inline const meta_type& currentMeta(void) const{ return(this->meta[this->metaPointer]); }
	inline const_reference operator[](const U32 p) const{ return this->runs[this->runsPointer + p]; }

	const U16 size(void) const{ return this->support->variants; }
	void reset(void){
		this->metaPointer = 0;
		this->runsPointer = 0;
	}

	void WriteVariant(const Totempole::TotempoleReader& totempole, IO::BasicBuffer& buffer, bool dropGenotypes = false) const{
		// All genotypes in this line will have the same phase
		const char separator = this->currentMeta().phased == 1 ? '|' : '/';

		// Note:
		// Much faster to first write to a char buffer then flush
		// instead of keep writing to cout (even without manual flushing)
		buffer += totempole.getContig(this->support->contigID).name;
		buffer += '\t';
		buffer += std::to_string(this->currentMeta().position);
		buffer += '\t';
		buffer += '.';
		buffer += '\t';
		buffer += Constants::REF_ALT_LOOKUP[this->currentMeta().ref_alt >> 4];
		buffer += '\t';
		buffer += Constants::REF_ALT_LOOKUP[this->currentMeta().ref_alt & ((1 << 4) - 1)];
		buffer += '\t';
		buffer += Constants::QUAL;
		buffer += '\t';
		buffer += Constants::PASS;
		buffer += '\t';
		buffer += std::string("HWE_P=");
		buffer += std::to_string(this->currentMeta().HWE_P);
		buffer += std::string(";MAF=");
		buffer += std::to_string(this->currentMeta().MAF);

		if(!dropGenotypes){
			buffer += '\t';
			buffer += Constants::GT;
			buffer += '\t';

			// For each run length encoded entry
			for(U32 i = 0; i < this->currentMeta().runs - 1; ++i){
				const char& left  = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[(*this)[i].alleleA];
				const char& right = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[(*this)[i].alleleB];

				// Repeat genotype run-length times
				for(U32 k = 0; k < (*this)[i].runs; ++k){
					buffer += left;
					buffer += separator;
					buffer += right;
					buffer += '\t';
				}
			}

			// For the last run length encoded entry
			const char& left  = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[(*this)[this->currentMeta().runs - 1].alleleA];
			const char& right = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[(*this)[this->currentMeta().runs - 1].alleleB];

			// Repeat genotype run-length - 1 times
			// Do not put a tab delimiter last
			for(U32 k = 0; k < (*this)[this->currentMeta().runs - 1].runs - 1; ++k){
				buffer += left;
				buffer += separator;
				buffer += right;
				buffer += '\t';
			}
			// Place a new line in the end instead
			buffer += left;
			buffer += separator;
			buffer += right;
			buffer += '\n';
		} else {
			buffer += '\n';
		}
	}

	bool buildPacked(const U64& samples);
	void clearPacked(void){ delete this->packed; }

public:
	U32 metaPointer;
	U32 runsPointer;
	const Totempole::TotempoleEntry* const support; // parent Totempole information
	const meta_type* const meta;
	const type* const runs;
	Base::GenotypeContainerBitvector* packed;
};

template <class T, class Y>
bool TomahawkBlock<T, Y>::buildPacked(const U64& samples){
	//return(this->packed->Build(*this, samples));
	return(true);
}


template <class T>
class TomahawkBlockManager{
	typedef TomahawkBlockManager<T>    self_type;
	typedef TomahawkBlock<const T, Support::TomahawkRun<T>>     controller_type;
	typedef Totempole::TotempoleEntry  totempole_entry_type;
	typedef Totempole::TotempoleReader totempole_reader_type;

public:
	TomahawkBlockManager(const totempole_reader_type& totempole_reader) :
		header(header)
	{}
	~TomahawkBlockManager(){}

	controller_type operator[](const U32 p) const{ return(controller_type(this->blocks[p])); } // copy constructor return
	void Add(const char* data, const totempole_entry_type& entry){ this->blocks.push_back(controller_type(data, entry)); }
	bool BuildVectorized(void){
		for(U32 i = 0; i < this->blocks.size(); ++i)
			this->blocks[i].buildPacked(header.getSamples());

		return true;
	}

	inline size_t size(void) const{ return this->blocks.size(); }
	U32 getVariants(void) const{
		U32 variants = 0;
		std::cerr << "test: " << this->size() << std::endl;

		for(U32 i = 0; i < this->size(); ++i){
			std::cerr << "debug: " << i << '\t' << this->size() << std::endl;
			std::cerr << &this->blocks[i] << std::endl;
			//std::cerr << this->blocks[i].
			variants += this->blocks[i].size();
		}

		std::cerr << "variants: " << variants << std::endl;

		return variants;
	}

public:
	std::vector<controller_type> blocks;

	// Header
	const Totempole::TotempoleReader& header;
};


}

#endif /* TOMAHAWK_TomahawkBlockManager_H_ */
