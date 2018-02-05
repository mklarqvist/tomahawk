#ifndef TOMAHAWK_BASE_GENOTYPE_CONTAINER_BITVECTOR_H_
#define TOMAHAWK_BASE_GENOTYPE_CONTAINER_BITVECTOR_H_

#include "TomahawkSupport.h"
#include "genotype_bitvector.h"
#include "genotype_container_runlength.h"

namespace Tomahawk{
namespace Base{

/**<
 * Meta container for bit-packed blocks of genotypes
 */
class GenotypeContainerBitvector{
private:
    typedef GenotypeContainerBitvector self_type;
    typedef Base::GenotypeBitvector<>  value_type;
    typedef value_type&                reference;
    typedef const value_type&          const_reference;
    typedef value_type*                pointer;
    typedef const value_type*          const_pointer;
    typedef std::ptrdiff_t             difference_type;
    typedef std::size_t                size_type;
    typedef Totempole::TotempoleEntry  support_type;

public:
	GenotypeContainerBitvector() : n_entries(0), n_capacity(0), __entries(nullptr){}
	~GenotypeContainerBitvector(){
		// Cleanup
		for(std::size_t i = 0; i < this->n_entries; ++i)
			((this->__entries + i)->~GenotypeBitvector)();

		::operator delete[](static_cast<void*>(this->__entries));
	}

	GenotypeContainerBitvector(const GenotypeContainerRunlength<BYTE>& genotype_container, const U64& n_samples) :
		n_entries(0),
		n_capacity(this->n_entries),
		__entries(nullptr)
	{
		this->Build<BYTE>(genotype_container, n_samples);
	}

	GenotypeContainerBitvector(const GenotypeContainerRunlength<U16>& genotype_container, const U64& n_samples) :
		n_entries(0),
		n_capacity(this->n_entries),
		__entries(nullptr)
	{
		this->Build<U16>(genotype_container, n_samples);
	}

	GenotypeContainerBitvector(const GenotypeContainerRunlength<U32>& genotype_container, const U64& n_samples) :
		n_entries(0),
		n_capacity(this->n_entries),
		__entries(nullptr)
	{
		this->Build<U32>(genotype_container, n_samples);
	}

	GenotypeContainerBitvector(const GenotypeContainerRunlength<U64>& genotype_container, const U64& n_samples) :
		n_entries(0),
		n_capacity(this->n_entries),
		__entries(nullptr)
	{
		this->Build<U64>(genotype_container, n_samples);
	}

	// copy constructor
	GenotypeContainerBitvector(const self_type& other) :
		n_entries(other.n_entries),
		n_capacity(other.n_capacity),
		__entries(other.__entries)
	{

	}

	// move constructor
	GenotypeContainerBitvector(self_type&& other) noexcept :
		n_entries(other.n_entries),
		n_capacity(other.n_capacity),
		__entries(other.__entries)
	{
		other.__entries = nullptr;
	}

	GenotypeContainerBitvector& operator=(self_type&& other) noexcept{
		 // prevent self-move
		if(this != &other){
			this->n_entries = other.n_entries;
			this->n_capacity = other.n_capacity;
			// swap
		}
		return *this;
	}

	class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const{ return *ptr_; }
		pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	class const_iterator{
	private:
		typedef const_iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		const_reference operator*() const{ return *ptr_; }
		const_pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	// Element access
	inline reference at(const size_type& position){ return(this->__entries[position]); }
	inline const_reference at(const size_type& position) const{ return(this->__entries[position]); }
	inline reference operator[](const size_type& position){ return(this->__entries[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->__entries[position]); }
	inline pointer data(void){ return(this->__entries); }
	inline const_pointer data(void) const{ return(this->__entries); }
	inline reference front(void){ return(this->__entries[0]); }
	inline const_reference front(void) const{ return(this->__entries[0]); }
	inline reference back(void){ return(this->__entries[this->n_entries - 1]); }
	inline const_reference back(void) const{ return(this->__entries[this->n_entries - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }
	inline const size_type& capacity(void) const{ return(this->n_capacity); }

	// Iterator
	inline iterator begin(){ return iterator(&this->__entries[0]); }
	inline iterator end()  { return iterator(&this->__entries[this->n_entries - 1]); }
	inline const_iterator begin()  const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->__entries[this->n_entries - 1]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->__entries[this->n_entries - 1]); }

	// Internal construction function
	template <class T>
	bool Build(const GenotypeContainerRunlength<T>& genotype_container, const U64& n_samples);

public:
	size_type n_entries;
	size_type n_capacity;
	pointer   __entries;
};

template <class T>
bool GenotypeContainerBitvector::Build(const GenotypeContainerRunlength<T>& genotype_container,
                                                                 const U64& n_samples)
{
	if(genotype_container.size() == 0)
		return false;

	// Cleanup
	for(std::size_t i = 0; i < this->n_entries; ++i)
		((this->__entries + i)->~GenotypeBitvector)();

	::operator delete[](static_cast<void*>(this->__entries));

	// Allocate new
	this->n_entries  = genotype_container.size();
	this->n_capacity = this->n_entries;
	this->__entries  = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	const U32 byte_width = ceil((double)n_samples/4);

	// INVERSE mask is cheaper in terms of instructions used
	// exploited in calculations: TomahawkCalculationSlave
	const BYTE lookup_mask[16] = {0, 0, 3, 3, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
	const BYTE lookup_data[16] = {0, 1, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	// Cycle over variants in container
	for(U32 i = 0; i < genotype_container.size(); ++i){
		new( &this->__entries[i] ) value_type( byte_width );
		Algorithm::GenotypeBitPacker packerA(this->__entries[i].data, 2);
		Algorithm::GenotypeBitPacker packerB(this->__entries[i].mask, 2);

		// Cycle over runs in container
		for(U32 j = 0; j < genotype_container[i].size(); ++j){
			const Support::TomahawkRunPacked<T>* const packed = reinterpret_cast<const Support::TomahawkRunPacked<T>* const>(&genotype_container[i][j]);
			packerA.add(lookup_data[packed->alleles], packed->runs);
			packerB.add(lookup_mask[packed->alleles], packed->runs);
		}
	}

	const U32 byteAlignedEnd  = byte_width / (GENOTYPE_TRIP_COUNT/4) * (GENOTYPE_TRIP_COUNT/4);

	// Search for zero runs in either end
	for(U32 i = 0; i < genotype_container.size(); ++i){
		S32 j = 0;

		// Search from left->right
		for(; j < byteAlignedEnd; ++j){
			if(this->__entries[i].data[j] != 0 || this->__entries[i].mask[j] != 0)
				break;
		}

		// Front of zeroes
		this->__entries[i].frontZero = ((j - 1 < 0 ? 0 : j - 1)*4)/GENOTYPE_TRIP_COUNT;
		if(j == byteAlignedEnd)
			continue;

		j = byteAlignedEnd - 1;
		for(; j > 0; --j){
			if(this->__entries[i].data[j] != 0 || this->__entries[i].mask[j] != 0)
				break;
		}

		// Tail of zeroes
		this->__entries[i].tailZero = ((byteAlignedEnd - (j+1))*4)/GENOTYPE_TRIP_COUNT;
	}
	return true;
}

}
}



#endif /* TOMAHAWK_BASE_GENOTYPE_CONTAINER_BITVECTOR_H_ */
