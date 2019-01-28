#include "index.h"

namespace tomahawk {

IndexEntry::IndexEntry() : rid(0), n(0), minpos(0), maxpos(0), b_unc(0), b_cmp(0), foff(0), fend(0){}
IndexEntry::~IndexEntry(){}

twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntry& self){
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

twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntry& self){
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

//
IndexEntryOutput::IndexEntryOutput() : ridB(-1){ rid = -1; }
IndexEntryOutput::~IndexEntryOutput(){}

void IndexEntryOutput::clear(){
	rid = -1; ridB = -1; minpos = 0; n = 0;
	foff = 0; fend = 0;
}

twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntryOutput& self){
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

twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntryOutput& self){
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

//
IndexEntryEntry::IndexEntryEntry() : nn(0){}

void IndexEntryEntry::operator+=(const IndexEntry& entry){
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

twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntryEntry& self){
	SerializePrimitive(self.rid, buffer);
	SerializePrimitive(self.n, buffer);
	SerializePrimitive(self.minpos, buffer);
	SerializePrimitive(self.maxpos, buffer);
	SerializePrimitive(self.foff, buffer);
	SerializePrimitive(self.fend, buffer);
	SerializePrimitive(self.nn, buffer);
	return(buffer);
}

twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntryEntry& self){
	DeserializePrimitive(self.rid, buffer);
	DeserializePrimitive(self.n, buffer);
	DeserializePrimitive(self.minpos, buffer);
	DeserializePrimitive(self.maxpos, buffer);
	DeserializePrimitive(self.foff, buffer);
	DeserializePrimitive(self.fend, buffer);
	DeserializePrimitive(self.nn, buffer);
	return(buffer);
}

// index
Index::Index(void) : n(0), m(0), m_ent(0), ent(nullptr), ent_meta(nullptr){}
Index::Index(const uint32_t n_contigs) : n(0), m(0), m_ent(n_contigs), ent(nullptr), ent_meta(new IndexEntryEntry[n_contigs]){}
Index::~Index(){ delete[] ent; delete[] ent_meta; }

void Index::Add(const IndexEntry& rec){
	if(this->n == this->m) this->resize(); // will also trigger when n=0 and m=0
	this->ent[this->n++] = rec;
	this->ent_meta[rec.rid] += rec;
}

//std::vector< IndexEntry* > FindOverlap(const uint32_t rid) const;
//std::vector< IndexEntry* > FindOverlap(const uint32_t rid, const uint32_t pos) const;
std::vector< IndexEntry* > Index::FindOverlap(const uint32_t rid, const uint32_t posA, const uint32_t posB) const{
	std::vector<IndexEntry*> ret;
	for(int i = 0; i < n; ++i){
		if(ent[i].rid == rid && ent[i].minpos <= posB && ent[i].maxpos >= posA){
			//std::cerr << "overlap=(" << rid << "," << posA << "," << posB << ") and (" << ent[i].rid << "," << ent[i].minpos << "," << ent[i].maxpos << ")" << std::endl;
			ret.push_back(&ent[i]);
		}
	}
	return(ret);
}


void Index::resize(void){
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
	m *= 2;
}

uint64_t Index::GetTotalVariants() const{
	uint64_t n_total = 0;
	for(int i = 0; i < n; ++i) n_total += ent[i].n;
	return(n_total);
}

twk_buffer_t& operator<<(twk_buffer_t& buffer, const Index& self){
	SerializePrimitive(TOMAHAWK_INDEX_START_MARKER, buffer);
	SerializePrimitive(self.n, buffer);
	SerializePrimitive(self.m, buffer);
	SerializePrimitive(self.m_ent, buffer);
	for(int i = 0; i < self.n; ++i) buffer << self.ent[i];
	for(int i = 0; i < self.m_ent; ++i) buffer << self.ent_meta[i];
	return(buffer);
}

twk_buffer_t& operator>>(twk_buffer_t& buffer, Index& self){
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

// inde output
IndexOutput::IndexOutput(void) : state(TWK_IDX_UNSORTED), n(0), m(0), m_ent(0), ent(nullptr), ent_meta(nullptr){}
IndexOutput::IndexOutput(const uint32_t n_contigs) : state(TWK_IDX_UNSORTED), n(0), m(0), m_ent(n_contigs), ent(nullptr), ent_meta(new IndexEntryEntry[n_contigs]){}
IndexOutput::~IndexOutput(){ delete[] ent; delete[] ent_meta; }

void IndexOutput::SetChroms(const uint32_t n_contigs){
	delete[] ent_meta;
	ent_meta = new IndexEntryEntry[n_contigs];
	m_ent = n_contigs;
}

void IndexOutput::Add(const IndexEntryOutput& rec){
	if(this->n == this->m) this->resize(); // will also trigger when n=0 and m=0
	this->ent[this->n++] = rec;
	//this->ent_meta[rec.rid] += rec;
}

void IndexOutput::AddThreadSafe(const IndexEntryOutput& rec){
	spinlock.lock();
	if(this->n == this->m) this->resize(); // will also trigger when n=0 and m=0
	this->ent[this->n++] = rec;
	//this->ent_meta[rec.rid] += rec;
	spinlock.unlock();
}

void IndexOutput::resize(void){
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

uint64_t IndexOutput::GetTotalVariants() const{
	uint64_t n_total = 0;
	for(int i = 0; i < n; ++i) n_total += ent[i].n;
	return(n_total);
}

//std::vector< IndexEntryOutput* > FindOverlap(const uint32_t rid) const;
//std::vector< IndexEntryOutput* > FindOverlap(const uint32_t rid, const uint32_t pos) const;
std::vector< IndexEntryOutput* > IndexOutput::FindOverlap(const uint32_t rid, const uint32_t posA, const uint32_t posB) const{
	std::vector<IndexEntryOutput*> ret;
	for(int i = 0; i < n; ++i){
		if(ent[i].rid == rid && ent[i].minpos <= posB && ent[i].maxpos >= posA){
			//std::cerr << "overlap=(" << rid << "," << posA << "," << posB << ") and (" << ent[i].rid << "," << ent[i].minpos << "," << ent[i].maxpos << ")" << std::endl;
			ret.push_back(&ent[i]);
		}
	}
	return(ret);
}

twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexOutput& self){
	SerializePrimitive(TOMAHAWK_INDEX_START_MARKER, buffer);
	SerializePrimitive(self.state, buffer);
	SerializePrimitive(self.n, buffer);
	SerializePrimitive(self.m, buffer);
	SerializePrimitive(self.m_ent, buffer);
	for(int i = 0; i < self.n; ++i) buffer << self.ent[i];
	for(int i = 0; i < self.m_ent; ++i) buffer << self.ent_meta[i];
	return(buffer);
}

twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexOutput& self){
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

}
