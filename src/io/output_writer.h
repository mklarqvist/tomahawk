#ifndef IO_OUTPUT_WRITER_H_
#define IO_OUTPUT_WRITER_H_

#include "../support/MagicConstants.h"
#include "../support/simd_definitions.h"
#include "../support/helpers.h"
#include "../io/compression/TGZFController.h"
#include "../algorithm/spinlock.h"
#include "../index/tomahawk_header.h"
#include "../index/index_entry.h"
#include "../tomahawk/output_container.h"
#include "../tomahawk/two/output_entry_support.h"
#include "../tomahawk/meta_entry.h"
#include "../index/index_entry.h"
#include "../index/index_container.h"
#include "../index/index.h"
#include "../index/footer.h"

namespace Tomahawk{
namespace IO{

class OutputWriter{
private:
	typedef OutputWriter                    self_type;
	typedef TGZFController                  compression_type;
	typedef Algorithm::SpinLock             spin_lock_type;
	typedef BasicBuffer                     buffer_type;
	typedef TomahawkHeader                  twk_header_type;
	typedef Totempole::IndexEntry           index_entry_type;
	typedef OutputEntry                     entry_type;
	typedef Support::OutputEntrySupport     entry_support_type;
	typedef Totempole::IndexEntry           header_entry_type;
	typedef Totempole::IndexContainer       index_container_type;
	typedef Index                           index_type;
	typedef Totempole::Footer               footer_type;
	typedef size_t                          size_type;
	typedef OutputContainer                 container_type;

public:
	OutputWriter(void);
	OutputWriter(std::string input_file);
	OutputWriter(const self_type& other);
	~OutputWriter(void);

	inline const U64& sizeEntries(void) const{ return(this->n_entries); }
	inline const U32& sizeBlocks(void) const{ return(this->n_blocks); }

	// Setters
	inline void setSorted(const bool yes){ this->writing_sorted_ = yes; }
	inline void setPartialSorted(const bool yes){ this->writing_sorted_partial_ = yes; }

	// Getters
	inline const bool isSorted(void) const{ return(this->writing_sorted_); }
	inline const bool isPartialSorted(void) const{ return(this->writing_sorted_partial_); }

	bool open(const std::string& output_file);
	int WriteHeaders(twk_header_type& twk_header);
	void WriteFinal(void);
	void flush(void);

	inline self_type& operator+=(const self_type& other){
		this->n_entries += other.n_entries;
		this->n_blocks  += other.n_blocks;
		if(other.l_largest_uncompressed > this->l_largest_uncompressed)
			this->l_largest_uncompressed = other.l_largest_uncompressed;

		return(*this);
	}

	inline self_type& operator=(const self_type& other){
		this->n_blocks  = other.n_blocks;
		this->n_entries = other.n_entries;
		if(other.l_largest_uncompressed > this->l_largest_uncompressed)
			this->l_largest_uncompressed = other.l_largest_uncompressed;
		return(*this);
	}

	template <class T>
	void Add(const MetaEntry<T>& meta_a, const MetaEntry<T>& meta_b, const header_entry_type& header_a, const header_entry_type& header_b, const entry_support_type& helper);

	inline void ResetProgress(void){ this->n_progress_count = 0; }
	inline const U32& getProgressCounts(void) const{ return this->n_progress_count; }

	inline void operator<<(const entry_type& entry){
		this->buffer << entry;
		++this->n_entries;

		if(this->buffer.size() > this->l_flush_limit)
			this->flush();
	}

	void operator<<(const container_type& container);
	void operator<<(buffer_type& buffer);

private:
	void CheckOutputNames(const std::string& input);

private:
	std::string      filename;
	std::string      basePath;
	std::string      baseName;
	bool             owns_pointers;
	bool             writing_sorted_;
	bool             writing_sorted_partial_;
	U64              n_entries;        // number of entries written
	U32              n_progress_count; // lines added since last flush
	U32              n_blocks;         // number of index blocks writtenflush_limit
	U32              l_flush_limit;
	U32              l_largest_uncompressed;
	index_entry_type index_entry;      // keep track of sort order
	std::ofstream*   stream;
	buffer_type      buffer;
	compression_type compressor;
	spin_lock_type*  spin_lock;
	index_type*      index_;
	footer_type*     footer_;
};

template <class T>
void OutputWriter::Add(const MetaEntry<T>& meta_a, const MetaEntry<T>& meta_b, const header_entry_type& header_a, const header_entry_type& header_b, const entry_support_type& helper){
	const U32 writePosA = meta_a.position << 2 | meta_a.phased << 1 | meta_a.missing;
	const U32 writePosB = meta_b.position << 2 | meta_b.phased << 1 | meta_b.missing;
	this->buffer += helper.controller;
	this->buffer += header_a.contigID;
	this->buffer += writePosA;
	this->buffer += header_b.contigID;
	this->buffer += writePosB;
	this->buffer << helper;
	// Add reverse
	this->buffer += helper.controller;
	this->buffer += header_b.contigID;
	this->buffer += writePosB;
	this->buffer += header_a.contigID;
	this->buffer += writePosA;
	this->buffer << helper;

	this->n_entries += 2;
	this->n_progress_count += 2;
	this->index_entry.n_variants += 2;

	if(this->buffer.size() > this->l_flush_limit)
		this->flush();
}

}
}

#endif /* IO_OUTPUT_WRITER_H_ */
