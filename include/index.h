/*
Copyright (C) 2016-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TWK_INDEX_H_
#define TWK_INDEX_H_

#include <unordered_map>
#include <cstdint>
#include <cassert>

#include "core.h"
#include "buffer.h"

namespace tomahawk {

struct IndexEntry {
	IndexEntry();
	virtual ~IndexEntry();

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntry& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntry& self);

public:
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
	IndexEntryOutput();
	virtual ~IndexEntryOutput();
	void clear();
	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntryOutput& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntryOutput& self);

public:
	int32_t ridB; // if ridB is mixed in this block we set this to -1
};

struct IndexEntryEntry : public IndexEntry {
	IndexEntryEntry();
	void operator+=(const IndexEntry& entry);
	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexEntryEntry& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexEntryEntry& self);

public:
	uint64_t nn; // number of consecutive blocks
};

class Index {
public:
	Index(void);
	Index(const uint32_t n_contigs);
	~Index();

	inline void operator+=(const IndexEntry& rec){ this->Add(rec); }

	void Add(const IndexEntry& rec);

	//std::vector< IndexEntry* > FindOverlap(const uint32_t rid) const;
	//std::vector< IndexEntry* > FindOverlap(const uint32_t rid, const uint32_t pos) const;
	std::vector< IndexEntry* > FindOverlap(const uint32_t rid, const uint32_t posA, const uint32_t posB) const;

	void resize(void);

	uint64_t GetTotalVariants() const;

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const Index& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, Index& self);

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
	IndexOutput(void);
	IndexOutput(const uint32_t n_contigs);
	~IndexOutput();

	void SetChroms(const uint32_t n_contigs);
	inline void operator+=(const IndexEntryOutput& rec){ this->Add(rec); }

	void Add(const IndexEntryOutput& rec);
	void AddThreadSafe(const IndexEntryOutput& rec);
	void resize(void);

	uint64_t GetTotalVariants() const;

	std::vector< IndexEntryOutput* > FindOverlap(const uint32_t rid, const uint32_t posA, const uint32_t posB) const;
	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const IndexOutput& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, IndexOutput& self);

public:
	uint8_t state; // sorted state of file
	uint64_t n, m, m_ent;
	IndexEntryOutput* ent;
	IndexEntryEntry* ent_meta;
	SpinLock spinlock;
};

}


#endif /* INDEX_H_ */
