#ifndef TOMAHAWK_BASE_GENOTYPE_CONTAINER_REFERENCE_H_
#define TOMAHAWK_BASE_GENOTYPE_CONTAINER_REFERENCE_H_

#include <cstring>  // size_t, ptrdiff_t

#include "../support/type_definitions.h"
#include "../index/index_entry.h"
#include "genotype_container_bitvector.h"
#include "haplotype_bitvector.h"

namespace Tomahawk{
namespace Base{

template <class T>
struct GenotypeRefEntry{
public:
	typedef GenotypeRefEntry<T>            self_type;
	typedef Support::GenotypeDiploidRun<T> genotype_type;
	typedef Base::GenotypeBitvector<>      genotype_bitvector_type;
	typedef MetaEntry                      meta_type;
	typedef Base::HaplotypeBitVector       haplotype_bitvector_type;

public:
	GenotypeRefEntry(const meta_type& meta_entry, const T* const genotypes) :
		meta_entry(meta_entry),
		genotypes(reinterpret_cast<genotype_type*>(genotypes)),
		genotype_bitvector(nullptr),
		haplotype_bitvector(nullptr)
	{

	}

	~GenotypeRefEntry(){
		delete this->genotype_bitvector;
		delete this->haplotype_bitvector;
	}

public:
	const meta_type& meta_entry;
	const genotype_type* const genotypes;
	genotype_bitvector_type* genotype_bitvector;
	haplotype_bitvector_type* haplotype_bitvector;
};

/**<
 * Special genotype container for both run-length encoded
 * and bit-encoded genotypes using unaligned memory directly
 * interpreted from type-casts at compile time and has some
 * psuedo-iterator capabilities.
 * This works because there is no random access in the pairwise
 * comparator functions. Upper triangular comparisons can be done
 * by strictly iterating forward.
 *
 * Data can only be interpreted by invoking the standard constructor.
 * All other constructors have non-standard meaning.
 * 1) Copy constructor: copies the pointer addresses and iterator
 *                      positions only.
 * 2) Assignment operator: copies the iterator position only
 */
template <class T>
class GenotypeContainerReference{
private:
	typedef GenotypeContainerReference  self_type;

protected:
	typedef Totempole::IndexEntry          header_entry_type;
	typedef GenotypeContainerBitvector     container_bitvector_type;
	typedef Base::GenotypeBitvector<>      genotype_bitvector_type;
	typedef Base::HaplotypeBitVector       haplotype_bitvector_type;
	typedef Support::GenotypeDiploidRun<T> value_type;
	typedef value_type&                    reference;
	typedef const value_type&              const_reference;
	typedef value_type*                    pointer;
	typedef const value_type*              const_pointer;
	typedef std::ptrdiff_t                 difference_type;
	typedef std::size_t                    size_type;
	typedef MetaEntry                      meta_type;

public:
	GenotypeContainerReference() :
		n_entries(0),
		iterator_position_meta(0),
		iterator_position_runs(0),
		build_bitvectors(false),
		owns_data(true),
		meta_entries(nullptr),
		genotype_entries(nullptr),
		index_entry(nullptr),
		bit_vectors(nullptr),
		haplotype_bitvectors(nullptr)
	{
	}

	GenotypeContainerReference(const char* const data,
	                           const size_t l_data,
	                           const header_entry_type& index_entry,
	                           const size_t n_samples,
	                           const bool build_bitvectors = true) :
		n_entries(index_entry.n_variants),
		iterator_position_meta(0),
		iterator_position_runs(0),
		build_bitvectors(build_bitvectors),
		owns_data(true),
		meta_entries(static_cast<meta_type*>(::operator new[](this->size()*sizeof(meta_type)))),
		genotype_entries(reinterpret_cast<const_pointer>(&data[this->size()*TOMAHAWK_ENTRY_META_SIZE])),
		index_entry(&index_entry),
		bit_vectors(nullptr),
		haplotype_bitvectors(nullptr)
	{
		if(l_data == 0) return;

		// Interpret meta entries
		size_t cumulative_position = 0;
		size_t genotype_cost = 0;

		for(size_t i = 0; i < this->size(); ++i){
			new( &this->meta_entries[i] ) meta_type( &data[cumulative_position] );
			cumulative_position += TOMAHAWK_ENTRY_META_SIZE;
			genotype_cost += meta_entries[i].runs*sizeof(T);
		}
		assert(cumulative_position + genotype_cost == l_data);

		if(build_bitvectors){
			this->bit_vectors = new container_bitvector_type();
			this->bit_vectors->Build(this->genotype_entries, this->meta_entries, this->size(), n_samples);

			this->haplotype_bitvectors = static_cast<haplotype_bitvector_type*>(::operator new[](this->size()*sizeof(haplotype_bitvector_type)));


			U32* tempList = new U32[n_samples*2]; // temporary vector
			U32 tempListIdx = 0;
			U32 cumulative_position = 0;
			for(U32 i = 0; i < this->size(); ++i){
				new( &this->haplotype_bitvectors[i] ) haplotype_bitvector_type( n_samples*2 );


				tempListIdx = 0;
				U32 cumsum = 0;
				for(U32 j = 0; j < meta_entries[i].runs; ++j, cumulative_position++){
					const Support::GenotypeDiploidRun<T>* const packed = reinterpret_cast<const Support::GenotypeDiploidRun<T>* const>(&this->genotype_entries[cumulative_position]);

					if((packed->alleleA & 3) != 0 || (packed->alleleB & 3) != 0){
						for(U32 k = 0; k < 2*packed->runs; k+=2){
							if((packed->alleleA & 3) != 0) tempList[tempListIdx++] = cumsum+k;
							if((packed->alleleB & 3) != 0) tempList[tempListIdx++] = cumsum+k+1;
							this->haplotype_bitvectors[i].l_list += (packed->alleleA & 3) != 0;
							this->haplotype_bitvectors[i].l_list += (packed->alleleB & 3) != 0;
							this->haplotype_bitvectors[i].set(cumsum+k,   (packed->alleleA & 3) != 0);
							this->haplotype_bitvectors[i].set(cumsum+k+1, (packed->alleleB & 3) != 0);
						}
					}
					cumsum += 2*packed->runs;
				}

				// Update index
				this->haplotype_bitvectors[i].indices = new U32[tempListIdx];
				for(U32 j = 0; j < tempListIdx; ++j)
					this->haplotype_bitvectors[i].indices[j] = tempList[j];
				//std::cerr << std::endl;
			}

			delete [] tempList;
		}
	}

	// Copy ctor: copies iterator positions and pointers
	GenotypeContainerReference(const self_type& other) :
		n_entries(other.n_entries),
		iterator_position_meta(other.iterator_position_meta),
		iterator_position_runs(other.iterator_position_runs),
		build_bitvectors(other.build_bitvectors),
		owns_data(false),
		meta_entries(other.meta_entries),
		genotype_entries(other.genotype_entries),
		index_entry(other.index_entry),
		bit_vectors(other.bit_vectors),
		haplotype_bitvectors(other.haplotype_bitvectors)
	{

	}

	~GenotypeContainerReference(){
		if(this->owns_data){
			if(this->build_bitvectors){
				delete this->bit_vectors;

				for(std::size_t i = 0; i < this->size(); ++i)
					(this->haplotype_bitvectors + i)->~HaplotypeBitVector();

				::operator delete[](static_cast<void*>(this->haplotype_bitvectors));
			}

			for(size_type i = 0; i < this->size(); ++i)
				((this->meta_entries + i)->~MetaEntry)();

			::operator delete[](static_cast<void*>(this->meta_entries));
		}
	}

	// // copy pointers only!
	void operator=(const self_type& other){
		this->iterator_position_runs = other.iterator_position_runs;
		this->iterator_position_meta = other.iterator_position_meta;
	}

	// Accessor
	inline const header_entry_type& getTotempole(void) const{ return(*this->index_entry); }
	inline const genotype_bitvector_type& getBitvector(const U32& position) const{ return(this->bit_vectors->at(position)); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }
	inline void resetIterator(void){ this->iterator_position_runs = 0; this->iterator_position_meta = 0; }

	// Accessor
	inline const char* const data(void) const{ return(this->genotype_entries); }
	inline const char* const meta_data(void) const{ return(this->meta_entries); }

	inline const meta_type& getMeta(const U32& position) const{ return(this->meta_entries[position]); }
	inline const meta_type& currentMeta(void) const{ return(this->meta_entries[this->iterator_position_meta]); }
	inline const meta_type& lastMeta(void) const{ return(this->meta_entries[this->n_entries - 1]); }
	inline const meta_type& firstMeta(void) const{ return(this->meta_entries[0]); }
	inline const_pointer current(void) const{ return(&this->genotype_entries[this->iterator_position_runs]); }
	inline const genotype_bitvector_type& currentBitvector(void) const{ return(this->bit_vectors->at(this->iterator_position_meta)); }
	inline const haplotype_bitvector_type& currentHaplotypeBitvector(void) const{ return(this->haplotype_bitvectors[this->iterator_position_meta]); }


	// Psuedo-iterator functionality
	inline void operator++(void){ this->iterator_position_runs += this->currentMeta().runs; ++this->iterator_position_meta; }
	inline void operator--(void){ this->iterator_position_runs -= this->currentMeta().runs; ++this->iterator_position_meta; }
	inline const_reference operator[](const U32& position) const{ return(this->genotype_entries[this->iterator_position_runs + position]); }
	inline const_reference at(const U32& position) const{ return(this->genotype_entries[this->iterator_position_runs + position]); }

protected:
	size_type                 n_entries;
	size_type                 iterator_position_meta;
	size_type                 iterator_position_runs;
	bool                      build_bitvectors;
	bool                      owns_data;
	meta_type*                meta_entries;
	const_pointer             genotype_entries;
	const header_entry_type*  index_entry;
	container_bitvector_type* bit_vectors;
	haplotype_bitvector_type* haplotype_bitvectors;
};

}
}

#endif /* TOMAHAWK_BASE_GENOTYPE_CONTAINER_REFERENCE_H_ */
