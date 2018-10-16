#ifndef TOMAHAWK_INDEX_H_
#define TOMAHAWK_INDEX_H_

#include <unordered_map>
#include <cstdint>
#include <cassert>

#include "core.h"
#include "buffer.h"

namespace tomahawk {

struct IndexEntry {
	IndexEntry() : rid(0), n(0), minpos(0), maxpos(0), foff(0), fend(0){}
	virtual ~IndexEntry(){}

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntry& self){
		SerializePrimitive(self.rid, buffer);
		SerializePrimitive(self.n, buffer);
		SerializePrimitive(self.minpos, buffer);
		SerializePrimitive(self.maxpos, buffer);
		SerializePrimitive(self.foff, buffer);
		SerializePrimitive(self.fend, buffer);
		return(buffer);
	}

	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntry& self){
		DeserializePrimitive(self.rid, buffer);
		DeserializePrimitive(self.n, buffer);
		DeserializePrimitive(self.minpos, buffer);
		DeserializePrimitive(self.maxpos, buffer);
		DeserializePrimitive(self.foff, buffer);
		DeserializePrimitive(self.fend, buffer);
		return(buffer);
	}

	uint32_t rid, n, minpos, maxpos;
	uint64_t foff, fend;
};

struct IndexEntryEntry : public IndexEntry {
	IndexEntryEntry() : nn(0){}

	void operator+=(const IndexEntry& entry){
		if(n == 0){
			minpos = entry.minpos;
			foff = entry.foff;
			rid = entry.rid;
		}
		assert(entry.rid == this->rid);

		n += entry.n;
		maxpos = entry.maxpos;
		fend = entry.fend;
		++nn;
	}

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntryEntry& self){
		SerializePrimitive(self.rid, buffer);
		SerializePrimitive(self.n, buffer);
		SerializePrimitive(self.minpos, buffer);
		SerializePrimitive(self.maxpos, buffer);
		SerializePrimitive(self.foff, buffer);
		SerializePrimitive(self.fend, buffer);
		SerializePrimitive(self.nn, buffer);
		return(buffer);
	}

	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntryEntry& self){
		DeserializePrimitive(self.rid, buffer);
		DeserializePrimitive(self.n, buffer);
		DeserializePrimitive(self.minpos, buffer);
		DeserializePrimitive(self.maxpos, buffer);
		DeserializePrimitive(self.foff, buffer);
		DeserializePrimitive(self.fend, buffer);
		DeserializePrimitive(self.nn, buffer);
		return(buffer);
	}

	uint32_t nn; // number of consecutive blocks
};

class Index {
public:
	Index(void) : n(0), m(0), m_ent(0), ent(nullptr), ent_meta(nullptr){}
	Index(const uint32_t n_contigs) : n(0), m(0), m_ent(n_contigs), ent(nullptr), ent_meta(new IndexEntryEntry[n_contigs]){}
	~Index(){ delete[] ent; delete[] ent_meta; }

	inline void operator+=(const IndexEntry& rec){ this->Add(rec); }

	void Add(const IndexEntry& rec){
		if(this->n == this->m) this->resize(); // will also trigger when n=0 and m=0
		this->ent[this->n++] = rec;
		this->ent_meta[rec.rid] += rec;
	}

	void resize(void){
		if(this->ent == nullptr){
			this->ent = new IndexEntry[500];
			this->n = 0;
			this->m = 500;
			return;
		}

		IndexEntry* temp = ent;
		ent = new IndexEntry[m*2];
		for(int i = 0; i < n; ++i) ent[i] = std::move(temp[i]); // move records over
		delete[] temp;
		m*=2;
	}

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const Index& self){
		SerializePrimitive(TOMAHAWK_INDEX_START_MARKER, buffer);
		SerializePrimitive(self.n, buffer);
		SerializePrimitive(self.m, buffer);
		SerializePrimitive(self.m_ent, buffer);
		for(int i = 0; i < self.n; ++i) buffer << self.ent[i];
		for(int i = 0; i < self.m_ent; ++i) buffer << self.ent_meta[i];
		return(buffer);
	}

	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, Index& self){
		uint64_t eof = 0;
		DeserializePrimitive(eof, buffer);
		assert(eof == TOMAHAWK_INDEX_START_MARKER);
		delete[] self.ent; delete[] self.ent_meta;
		DeserializePrimitive(self.n, buffer);
		DeserializePrimitive(self.m, buffer);
		DeserializePrimitive(self.m_ent, buffer);
		self.ent = new IndexEntry[self.m];
		self.ent_meta = new IndexEntryEntry[self.m_ent];
		for(int i = 0; i < self.n; ++i) buffer >> self.ent[i];
		for(int i = 0; i < self.m_ent; ++i) buffer >> self.ent_meta[i];
		return(buffer);
	}


public:
	uint32_t n, m, m_ent;
	IndexEntry* ent;
	IndexEntryEntry* ent_meta;
};

}


#endif /* INDEX_H_ */
