#ifndef TOMAHAWK_INDEX_H_
#define TOMAHAWK_INDEX_H_

#include <unordered_map>
#include <cstdint>
#include <cassert>

#include "core.h"
#include "buffer.h"
#include "spinlock.h"

namespace tomahawk {

struct IndexEntry {
	IndexEntry() : rid(0), n(0), minpos(0), maxpos(0), b_unc(0), b_cmp(0), foff(0), fend(0){}
	virtual ~IndexEntry(){}

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntry& self){
		SerializePrimitive(self.rid, buffer);
		SerializePrimitive(self.n, buffer);
		SerializePrimitive(self.minpos, buffer);
		SerializePrimitive(self.maxpos, buffer);
		SerializePrimitive(self.b_unc, buffer);
		SerializePrimitive(self.b_cmp, buffer);
		SerializePrimitive(self.foff, buffer);
		SerializePrimitive(self.fend, buffer);
		return(buffer);
	}

	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntry& self){
		DeserializePrimitive(self.rid, buffer);
		DeserializePrimitive(self.n, buffer);
		DeserializePrimitive(self.minpos, buffer);
		DeserializePrimitive(self.maxpos, buffer);
		DeserializePrimitive(self.b_unc, buffer);
		DeserializePrimitive(self.b_cmp, buffer);
		DeserializePrimitive(self.foff, buffer);
		DeserializePrimitive(self.fend, buffer);
		return(buffer);
	}

	int32_t rid;
	uint32_t n, minpos, maxpos;
	uint32_t b_unc, b_cmp;
	uint64_t foff, fend;
};

/**<
 * Index entry for two files. This entry is different as we require knowledge
 * of both the from rid:pos and the to rid as a tuple (rid:pos, rid).
 * This is different from the case in twk files where the tuple is limited to
 * (rid:pos) only.
 */
struct IndexEntryOutput : public IndexEntry {
	IndexEntryOutput() : ridB(-1){ rid = -1; }
	virtual ~IndexEntryOutput(){}

	void clear(){
		rid = -1; ridB = -1; minpos = 0; n = 0;
		foff = 0; fend = 0;
	}

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntryOutput& self){
		SerializePrimitive(self.rid, buffer);
		SerializePrimitive(self.n, buffer);
		SerializePrimitive(self.minpos, buffer);
		SerializePrimitive(self.maxpos, buffer);
		SerializePrimitive(self.b_unc, buffer);
		SerializePrimitive(self.b_cmp, buffer);
		SerializePrimitive(self.foff, buffer);
		SerializePrimitive(self.fend, buffer);
		SerializePrimitive(self.ridB, buffer);
		return(buffer);
	}

	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntryOutput& self){
		DeserializePrimitive(self.rid, buffer);
		DeserializePrimitive(self.n, buffer);
		DeserializePrimitive(self.minpos, buffer);
		DeserializePrimitive(self.maxpos, buffer);
		DeserializePrimitive(self.b_unc, buffer);
		DeserializePrimitive(self.b_cmp, buffer);
		DeserializePrimitive(self.foff, buffer);
		DeserializePrimitive(self.fend, buffer);
		DeserializePrimitive(self.ridB, buffer);
		return(buffer);
	}

public:
	int32_t ridB; // if ridB is mixed in this block we set this to -1
};

struct IndexEntryEntry : public IndexEntry {
	IndexEntryEntry() : nn(0){}

	void operator+=(const IndexEntry& entry){
		if(n == 0){
			minpos = entry.minpos;
			foff = entry.foff;
			rid = entry.rid;
		}

		// Assertion.
		if(entry.rid != this->rid){
			std::cerr << "illegal addition of rid: " << entry.rid << "!=" << rid << std::endl;
			std::cerr << entry.rid << ":" << entry.minpos << "-" << entry.maxpos << " " << entry.n << std::endl;
			exit(1);
		}

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

public:
	uint64_t nn; // number of consecutive blocks
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

	std::vector< IndexEntry* > FindOverlap(const uint32_t rid) const;
	std::vector< IndexEntry* > FindOverlap(const uint32_t rid, const uint32_t pos) const;
	std::vector< IndexEntry* > FindOverlap(const uint32_t rid, const uint32_t posA, const uint32_t posB) const{
		std::vector<IndexEntry*> ret;
		for(int i = 0; i < n; ++i){
			if(ent[i].rid == rid && ent[i].minpos <= posB && ent[i].maxpos >= posA){
				//std::cerr << "overlap=(" << rid << "," << posA << "," << posB << ") and (" << ent[i].rid << "," << ent[i].minpos << "," << ent[i].maxpos << ")" << std::endl;
				ret.push_back(&ent[i]);
			}
		}
		return(ret);
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

	uint64_t GetTotalVariants() const{
		uint64_t n_total = 0;
		for(int i = 0; i < n; ++i) n_total += ent[i].n;
		return(n_total);
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
	uint64_t n, m, m_ent;
	IndexEntry* ent;
	IndexEntryEntry* ent_meta;
};

#define TWK_IDX_UNSORTED 0
#define TWK_IDX_PARTIAL  1
#define TWK_IDX_SORTED   2

class IndexOutput {
public:
	IndexOutput(void) : state(TWK_IDX_UNSORTED), n(0), m(0), m_ent(0), ent(nullptr), ent_meta(nullptr){}
	IndexOutput(const uint32_t n_contigs) : state(TWK_IDX_UNSORTED), n(0), m(0), m_ent(n_contigs), ent(nullptr), ent_meta(new IndexEntryEntry[n_contigs]){}
	~IndexOutput(){ delete[] ent; delete[] ent_meta; }

	void SetChroms(const uint32_t n_contigs){
		delete[] ent_meta;
		ent_meta = new IndexEntryEntry[n_contigs];
		m_ent = n_contigs;
	}
	inline void operator+=(const IndexEntryOutput& rec){ this->Add(rec); }

	void Add(const IndexEntryOutput& rec){
		if(this->n == this->m) this->resize(); // will also trigger when n=0 and m=0
		this->ent[this->n++] = rec;
		//this->ent_meta[rec.rid] += rec;
	}

	void AddThreadSafe(const IndexEntryOutput& rec){
		spinlock.lock();
		if(this->n == this->m) this->resize(); // will also trigger when n=0 and m=0
		this->ent[this->n++] = rec;
		//this->ent_meta[rec.rid] += rec;
		spinlock.unlock();
	}

	void resize(void){
		if(this->ent == nullptr){
			this->ent = new IndexEntryOutput[500];
			this->n = 0;
			this->m = 500;
			return;
		}

		IndexEntryOutput* temp = ent;
		ent = new IndexEntryOutput[m*2];
		for(int i = 0; i < n; ++i) ent[i] = std::move(temp[i]); // move records over
		delete[] temp;
		m*=2;
	}

	uint64_t GetTotalVariants() const{
		uint64_t n_total = 0;
		for(int i = 0; i < n; ++i) n_total += ent[i].n;
		return(n_total);
	}

	//std::vector< IndexEntryOutput* > FindOverlap(const uint32_t rid) const;
	//std::vector< IndexEntryOutput* > FindOverlap(const uint32_t rid, const uint32_t pos) const;
	std::vector< IndexEntryOutput* > FindOverlap(const uint32_t rid, const uint32_t posA, const uint32_t posB) const{
		std::vector<IndexEntryOutput*> ret;
		for(int i = 0; i < n; ++i){
			if(ent[i].rid == rid && ent[i].minpos <= posB && ent[i].maxpos >= posA){
				//std::cerr << "overlap=(" << rid << "," << posA << "," << posB << ") and (" << ent[i].rid << "," << ent[i].minpos << "," << ent[i].maxpos << ")" << std::endl;
				ret.push_back(&ent[i]);
			}
		}
		return(ret);
	}

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexOutput& self){
		SerializePrimitive(TOMAHAWK_INDEX_START_MARKER, buffer);
		SerializePrimitive(self.state, buffer);
		SerializePrimitive(self.n, buffer);
		SerializePrimitive(self.m, buffer);
		SerializePrimitive(self.m_ent, buffer);
		for(int i = 0; i < self.n; ++i) buffer << self.ent[i];
		for(int i = 0; i < self.m_ent; ++i) buffer << self.ent_meta[i];
		return(buffer);
	}

	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexOutput& self){
		uint64_t eof = 0;
		DeserializePrimitive(eof, buffer);
		assert(eof == TOMAHAWK_INDEX_START_MARKER);
		delete[] self.ent; delete[] self.ent_meta;
		DeserializePrimitive(self.state, buffer);
		DeserializePrimitive(self.n, buffer);
		DeserializePrimitive(self.m, buffer);
		DeserializePrimitive(self.m_ent, buffer);
		self.ent = new IndexEntryOutput[self.m];
		self.ent_meta = new IndexEntryEntry[self.m_ent];
		for(int i = 0; i < self.n; ++i) buffer >> self.ent[i];
		for(int i = 0; i < self.m_ent; ++i) buffer >> self.ent_meta[i];
		return(buffer);
	}

public:
	uint8_t state; // sorted state of file
	uint64_t n, m, m_ent;
	IndexEntryOutput* ent;
	IndexEntryEntry* ent_meta;
	SpinLock spinlock;
};

}


#endif /* INDEX_H_ */
