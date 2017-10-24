#ifndef TOMAHAWK_TomahawkBlockManager_H_
#define TOMAHAWK_TomahawkBlockManager_H_

#include "../support/simd_definitions.h"
#include "../algorithm/GenotypeBitPacker.h"
#include "../tomahawk/base/TomahawkEntryMeta.h"
#include "../totempole/TotempoleEntry.h"
#include "../totempole/TotempoleReader.h"
#include "base/TomahawkRun.h"

namespace Tomahawk{

template <class T, class Y = Support::TomahawkRun<T>>
struct TomahawkIterator; // forward declare: required for build function

template <int T = SIMD_ALIGNMENT>
struct TomahawkBlockPackedPair{
public:
	TomahawkBlockPackedPair(const U32 size):
		frontZero(0),
		tailZero(0),
		frontZeroMissing(0),
		tailZeroMissing(0),
	#if SIMD_AVAILABLE == 1
		data((BYTE*)_mm_malloc(size, T)),
		mask((BYTE*)_mm_malloc(size, T))
	#else
		data(new BYTE[size]),
		mask(new BYTE[size])
	#endif
	{
		memset(this->data, 0, size);
		memset(this->mask, 0, size);
	}

	~TomahawkBlockPackedPair(){
	#if SIMD_AVAILABLE == 1
		_mm_free(this->data);
		_mm_free(this->mask);
	#else
		delete [] this->data;
		delete [] this->mask;
	#endif
	}

public:
	U32 frontZero;				// leading zeros in aligned vector width
	U32 tailZero;				// trailing zeros in aligned vector width
	U32 frontZeroMissing;		// number of missing values in leading zeros
	U32 tailZeroMissing;		// number of missing values in trailing zeros
	BYTE* data;
	BYTE* mask;
} __attribute__((aligned(16)));

class TomahawkBlockPacked{
	typedef TomahawkBlockPackedPair<> pair_type;

public:
	TomahawkBlockPacked() : width(0), data(nullptr){}
	~TomahawkBlockPacked(){
		delete [] this->data;
		delete this->data;
	}

	// copy constructor
	TomahawkBlockPacked(const TomahawkBlockPacked& other) :
		width(other.width),
		data(other.data)
	{

	}

	// move constructor
	TomahawkBlockPacked(TomahawkBlockPacked&& other) noexcept :
		width(other.width),
		data(other.data)
	{
		other.data = nullptr;
	}

	/** Move assignment operator */
	TomahawkBlockPacked& operator=(TomahawkBlockPacked&& other) noexcept{
		 // prevent self-move
		if(this != &other)
			this->width = other.width;

		return *this;
	}

	template <class T>
	bool Build(TomahawkIterator<T>& controller, const U64& samples);

	inline const pair_type& getData(const U32 p) const{ return(*this->data[p]); }

public:
	U32 width;
	pair_type** data;
};


template <class T, class Y>
struct TomahawkIterator{
	typedef Support::TomahawkEntryMetaRLE<T> meta_type;
	typedef Totempole::TotempoleEntry totempole_type;

public:
	TomahawkIterator(const char* target, const totempole_type& support) :
		position(0),
		metaPointer(0),
		runsPointer(0),
		support(&support),
		meta(reinterpret_cast<const meta_type* const>(target)),
		runs(reinterpret_cast<const Y* const>(&target[(TOMAHAWK_ENTRY_META_SIZE + sizeof(T)) * support.variants])),
		packed(new TomahawkBlockPacked)
	{

	}

	~TomahawkIterator() noexcept{}

	// copy constructor
	TomahawkIterator(const TomahawkIterator& other) :
		position(other.position),
		metaPointer(other.metaPointer),
		runsPointer(other.runsPointer),
		support(other.support),
		meta(other.meta),
		runs(other.runs),
		packed(other.packed)
	{

	}

	void operator=(const TomahawkIterator& other){
		this->position = other.position;
		this->metaPointer = other.metaPointer;
		this->runsPointer = other.runsPointer;
	}

	// move constructor
	TomahawkIterator(TomahawkIterator&& other) noexcept :
		position(other.position),
		metaPointer(other.metaPointer),
		runsPointer(other.runsPointer),
		support(other.support),
		meta(other.meta),
		runs(other.runs),
		packed(other.packed)
	{

	}

	inline void updatePacked(const TomahawkIterator& self){
		this->packed = new TomahawkBlockPacked(self.packed);
	}


	inline void operator++(void){
		++this->position;
		this->runsPointer += this->meta[this->metaPointer].runs;
		++this->metaPointer;
	}

	inline void operator--(void){
		--this->position;
		--this->metaPointer;
		this->runsPointer -= this->meta[this->metaPointer].runs;
	}

	bool nextVariant(const Support::TomahawkRun<T>*& run, const meta_type*& meta){
		if(this->position == this->support->variants)
			return false;

		run = &runs[this->runsPointer];
		meta = &this->meta[this->metaPointer];
		++(*this);

		return true;
	}

	bool nextVariant(const Support::TomahawkRunPacked<T>*& run, const meta_type*& meta){
		if(this->position == this->support->variants)
			return false;

		run = reinterpret_cast<const Support::TomahawkRunPacked<T>*>(&runs[this->runsPointer]);
		meta = &this->meta[this->metaPointer];
		++(*this);

		return true;
	}

	inline const meta_type& currentMeta(void) const{ return(this->meta[this->metaPointer]); }
	inline const Y& operator[](const U32 p) const{ return this->runs[this->runsPointer + p]; }

	const U16& size(void) const{ return this->support->variants; }
	void reset(void){
		this->metaPointer = 0;
		this->runsPointer = 0;
	}

	void WriteVariant(const Totempole::TotempoleReader& totempole, IO::BasicBuffer& buffer, bool dropGenotypes = false) const{
		// All genotypes in this line will have the same phase
		//const char separator = this->currentMeta().phased == 1 ? '|' : '/';

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
		buffer += std::string(";MGF=");
		buffer += std::to_string(this->currentMeta().MGF);

		if(!dropGenotypes){
			buffer += '\t';
			buffer += Constants::GT;
			buffer += '\t';

			// For each run length encoded entry
			for(U32 i = 0; i < this->currentMeta().runs - 1; ++i){
				const char& left  = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[(*this)[i].alleleA];
				const char& right = Constants::TOMAHAWK_ALLELE_LOOKUP_REVERSE[(*this)[i].alleleB];

				// Repeat genotype run-length times
				const char separator = (*this)[i].phasing == 1 ? '|' : '/';
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
			const char separator = (*this)[this->currentMeta().runs - 1].phasing == 1 ? '|' : '/';
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
	U32 position;
	U32 metaPointer;
	U32 runsPointer;
	const totempole_type* const support; // parent Totempole information
	const meta_type* const meta;
	const Y* const runs;
	TomahawkBlockPacked* packed;
};

template <class T>
bool TomahawkBlockPacked::Build(TomahawkIterator<T>& controller, const U64& samples){
	if(controller.support->variants == 0)
		return false;

	controller.reset();
	TomahawkIterator<T, Support::TomahawkRunPacked<T> >& c = *reinterpret_cast<TomahawkIterator<T, Support::TomahawkRunPacked<T> >*>(&controller);

	this->width = c.support->variants;
	this->data = new pair_type*[c.support->variants]; // Todo: fix this, poor locality

	const U32 byte_width = ceil((double)samples/4);

	// INVERSE mask is cheaper in terms of instructions used
	// exploited in calculations: TomahawkCalculationSlave
	const BYTE lookup_mask[16] = {0, 0, 3, 3, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
	const BYTE lookup_data[16] = {0, 1, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	for(U32 i = 0; i < c.support->variants; ++i){
		this->data[i] = new pair_type(byte_width);
		Algorithm::GenotypeBitPacker packerA(this->data[i]->data, 2);
		Algorithm::GenotypeBitPacker packerB(this->data[i]->mask, 2);

		for(U32 j = 0; j < c.meta[i].runs; ++j){
			packerA.add(lookup_data[c[j].alleles], c[j].runs);
			packerB.add(lookup_mask[c[j].alleles], c[j].runs);
		}
		++c;
	}
	controller.reset();

	const U32 byteAlignedEnd  = byte_width/(GENOTYPE_TRIP_COUNT/4)*(GENOTYPE_TRIP_COUNT/4);

	// Search for zero runs in either end
	for(U32 i = 0; i < c.support->variants; ++i){
		S32 j = 0;

		// Search from left->right
		for(; j < byteAlignedEnd; ++j){
			if(this->data[i]->data[j] != 0 || this->data[i]->mask[j] != 0)
				break;
		}

		// Front of zeroes
		this->data[i]->frontZero = ((j - 1 < 0 ? 0 : j - 1)*4)/GENOTYPE_TRIP_COUNT;
		if(j == byteAlignedEnd)
			break;

		j = byteAlignedEnd - 1;
		for(; j > 0; --j){
			if(this->data[i]->data[j] != 0 || this->data[i]->mask[j] != 0)
				break;
		}

		// Tail of zeroes
		this->data[i]->tailZero = ((byteAlignedEnd - (j+1))*4)/GENOTYPE_TRIP_COUNT;
	}
	return true;
}

template <class T, class Y>
bool TomahawkIterator<T, Y>::buildPacked(const U64& samples){
	return(this->packed->Build(*this, samples));
}


template <class T>
class TomahawkBlockManager{
	typedef TomahawkBlockManager<T> self_type;
	typedef TomahawkIterator<const T> controller_type;
	typedef Totempole::TotempoleEntry totempole_entry_type;

public:
	TomahawkBlockManager(const Totempole::TotempoleReader& header) :
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

		for(U32 i = 0; i < this->size(); ++i)
			variants += this->blocks[i].size();

		return variants;
	}

public:
	std::vector<controller_type> blocks;

	// Header
	const Totempole::TotempoleReader& header;
};


}

#endif /* TOMAHAWK_TomahawkBlockManager_H_ */
