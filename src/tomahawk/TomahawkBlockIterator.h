#ifndef TOMAHAWK_TOMAHAWKBLOCKITERATOR_H_
#define TOMAHAWK_TOMAHAWKBLOCKITERATOR_H_

#include "base/TomahawkEntryMeta.h"
#include "base/TomahawkGTEntries.h"

namespace Tomahawk{

template <class T>
class TomahawkBlockIterator{
	typedef TomahawkBlockIterator self_type;
	typedef Support::TomahawkEntryMeta<T> meta_type;
	typedef Support::TomahawkSupport meta_complex_type;
	typedef Totempole::TotempoleEntry totempole_entry_type;

public:
	TomahawkBlockIterator(char* data, const U64& size, const totempole_entry_type& totempole) :
		position(0),
		p_rle(0),
		p_simple(0),
		pointer(0),
		upper_limit(0),
		width(size),
		totempole(totempole),
		data(data),
		meta(reinterpret_cast<const meta_type*>(data)),
		encoding_RLE(&data[totempole.l_meta]),
		encoding_simple(&data[totempole.l_meta + totempole.l_gt_rle]),
		meta_complex(&data[totempole.l_meta + totempole.l_gt_rle + totempole.l_gt_simple])
	{
		this->upper_limit = this->meta[0].n_runs;
	}

	~TomahawkBlockIterator(){}

	bool operator++(void){
		if(this->position + 1 == this->totempole.n_variants) return false;
		++this->position;
		this->pointer = 0;

		this->upper_limit = this->getMeta().n_runs;
		if(this->getMeta().controller.biallelic){
			this->encoding_RLE = &this->data[totempole.l_meta + this->getMeta().virtual_offset_gt];
			++this->p_rle;
		}
		else {
			this->encoding_simple = &this->data[totempole.l_meta + totempole.l_gt_rle + this->getMeta().virtual_offset_gt];
			++this->p_simple;
		}

		return true;
	}

	inline const bool isRLE(void) const{ return(this->meta[this->position].isRLE()); }
	inline const U32& size(void) const{ return(this->upper_limit); }

	template <class S>
	const bool nextRun(const S*& run){
		if(this->pointer == this->upper_limit)
			return false;

		//run = this->encoding_RLE;
		run = reinterpret_cast<const S*>(this->encoding_RLE);
		++this->pointer;
		this->encoding_RLE += sizeof(S);
		return true;
	}

	template <class S, BYTE missing = 1>
	const bool nextRun(const Support::TomahawkRun<S, missing>*& run){
		if(this->pointer == this->upper_limit)
			return false;

		//run = this->encoding_RLE;
		run = reinterpret_cast<const Support::TomahawkRun<S, missing>*>(this->encoding_RLE);
		++this->pointer;
		this->encoding_RLE += sizeof(S);
		return true;
	}

	template <class S, BYTE missing = 1>
	const bool nextRun(const Support::TomahawkRunNoPhase<S, missing>*& run){
		if(this->pointer == this->upper_limit)
			return false;

		//run = this->encoding_RLE;
		run = reinterpret_cast<const Support::TomahawkRunNoPhase<S, missing>*>(this->encoding_RLE);
		++this->pointer;
		this->encoding_RLE += sizeof(S);
		return true;
	}

	template <class S>
	const bool nextRunSimple(const Support::TomahawkRunSimple<S>*& field){
		if(this->pointer == this->upper_limit)
			return false;

		field = reinterpret_cast<const Support::TomahawkRunSimple<S>*>(this->encoding_simple);
		++this->pointer;
		this->encoding_simple += sizeof(S);
		return true;
	}

	template <class S>
	const bool nextRunSimple(const S*& field){
		if(this->pointer == this->upper_limit)
			return false;

		field = reinterpret_cast<const S*>(this->encoding_simple);
		++this->pointer;
		this->encoding_simple += sizeof(S);
		return true;
	}

	inline const meta_type& getMeta(void) const{ return(this->meta[this->position]); }
	inline meta_complex_type& getMetaComplex(void){
		return(*reinterpret_cast<meta_complex_type*>(&this->meta_complex[this->meta[this->position].virtual_offset_cold_meta]));
	}

	bool countGenotypes(void);
	bool countGenotypesGroup(void);

private:
	bool __countGenotypesRLE(void);
	bool __countGenotypesRLEGroup(void);
	bool __countGenotypesSimple(void);
	bool __countGenotypesSimpleGroup(void);

private:
	U32 position;    // current meta position
	U32 p_rle;       // position RLE
	U32 p_simple;    // position simple
	U32 pointer;     // iterator pointer
	U32 upper_limit; // relative upper bounds in iterator
	const U64 width;
	const totempole_entry_type& totempole;
	const char* const data;
	const meta_type* meta;
	const char* encoding_RLE;
	const char* encoding_simple;
	char* meta_complex;
};

}

#endif /* TOMAHAWK_TOMAHAWKBLOCKITERATOR_H_ */
