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

/**<
 * Writer class for `two` entries. This class supports parallel writing
 * with the use of a lock-free spin-lock (requires C++11 because of the use
 * of atomic values). In parallel computing, each slave constructs their own
 * OutputWriter by invoking the copy-ctor and borrowing pointers from the
 * main instance.
 */
class OutputWriter{
protected:
	typedef OutputWriter                self_type;
	typedef TGZFController              compression_type;
	typedef Algorithm::SpinLock         spin_lock_type;
	typedef BasicBuffer                 buffer_type;
	typedef TomahawkHeader              twk_header_type;
	typedef Totempole::IndexEntry       index_entry_type;
	typedef OutputEntry                 entry_type;
	typedef Support::OutputEntrySupport entry_support_type;
	typedef Totempole::IndexEntry       header_entry_type;
	typedef Totempole::IndexContainer   index_container_type;
	typedef Index                       index_type;
	typedef Totempole::Footer           footer_type;
	typedef size_t                      size_type;
	typedef OutputContainer             container_type;

public:
	OutputWriter(void);
	OutputWriter(const self_type& other);
	virtual ~OutputWriter(void);

	inline const U64& sizeEntries(void) const{ return(this->n_entries); }
	inline const U32& sizeBlocks(void) const{ return(this->n_blocks); }

	// Setters
	inline void setSorted(const bool yes){ this->writing_sorted_ = yes; }
	inline void setPartialSorted(const bool yes){ this->writing_sorted_partial_ = yes; }
	inline void setFlushLimit(const U32 limit){ this->l_flush_limit = limit; }

	// Getters
	inline const bool isSorted(void) const{ return(this->writing_sorted_); }
	inline const bool isPartialSorted(void) const{ return(this->writing_sorted_partial_); }
	inline index_type& getIndex(void) const{ return(*this->index_); }
	inline entry_type& getPreviousEntry(void){ return(this->previous_entry); }
	inline index_entry_type& getCurrentIndexEntry(void){ return(this->index_entry); }

	virtual bool open(const std::string& output_file) =0;
	virtual int writeHeaders(twk_header_type& twk_header) =0;
	virtual void writeFinal(void) =0;
	virtual void flush(void) =0;

	inline void ResetProgress(void){ this->n_progress_count = 0; }
	inline const U32& getProgressCounts(void) const{ return this->n_progress_count; }

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

	/**<
	 * Primary function writing `two` entries to disk after being computed by a
	 * slave.
	 * @param meta_a   Meta information for the from container
	 * @param meta_b   Meta information for the to container
	 * @param header_a Tomahawk index entry for the from container
	 * @param header_b Tomahawk index entry for the to container
	 * @param helper   Helper structure used in computing LD. Holds the allele/genotype counts and statistics
	 */
	template <class T>
	void Add(const MetaEntry<T>& meta_a, const MetaEntry<T>& meta_b, const header_entry_type& header_a, const header_entry_type& header_b, const entry_support_type& helper);

	/**<
	 * Overloaded operator for adding a single `two` entry
	 * @param entry Input `two` entry
	 */
	inline void operator<<(const entry_type& entry){
		this->buffer << entry;
		++this->n_entries;

		// Check if the buffer has to be flushed after adding this entry
		if(this->buffer.size() > this->l_flush_limit)
			this->flush();
	}

	/**<
	 *
	 * @param entry
	 */
	inline void addSorted(const entry_type& entry){
		if(this->previous_entry.AcontigID != entry.AcontigID
			|| this->previous_entry.BcontigID != entry.BcontigID){
			//std::cerr << "switch in contig" << std::endl;
			this->flush();
			this->index_entry = entry;
		}

		// Check if the buffer has to be flushed after adding this entry
		if(this->buffer.size() > this->l_flush_limit){
			this->flush();
			this->index_entry = entry;
		}

		this->buffer << entry;
		++this->n_entries;
		this->previous_entry = entry;
		this->index_entry.max_position = entry.Aposition;
	}

	/**<
	 * Overloaded operator for adding an entire container of `two` entries
	 * @param container Target container of entries
	 */
	void operator<<(const container_type& container);

	/**<
	 * Overloaded operator for adding an entire buffer of `two` entries
	 * @param buffer Target buffer of entries
	 */
	virtual void operator<<(buffer_type& buffer) =0;

	virtual void writePrecompressedBlock(buffer_type& buffer, const U64& uncompressed_size) =0;

protected:
	virtual void CheckOutputNames(const std::string& input) =0;

protected:
	bool             owns_pointers;
	bool             writing_sorted_;
	bool             writing_sorted_partial_;
	U64              n_entries;        // number of entries written
	U32              n_progress_count; // lines added since last flush
	U32              n_blocks;         // number of index blocks writtenflush_limit
	U32              l_flush_limit;
	U32              l_largest_uncompressed;
	entry_type       previous_entry;
	index_entry_type index_entry;      // keep track of sort order

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

class OutputWriterStream : public OutputWriter{
private:
	typedef OutputWriterStream self_type;
	typedef OutputWriter       parent_type;

public:
	OutputWriterStream(void);
	OutputWriterStream(const self_type& other);
	~OutputWriterStream(void);

	bool open(const std::string& output_file){ return true; }
	int writeHeaders(twk_header_type& twk_header);
	void writeFinal(void);
	void flush(void);
	void operator<<(buffer_type& buffer);
	void operator<<(const container_type& container){ parent_type::operator<<(container); }
	void writePrecompressedBlock(buffer_type& buffer, const U64& uncompressed_size);

private:
	U64 b_written;
};

class OutputWriterFile : public OutputWriter {
private:
	typedef OutputWriterFile self_type;
	typedef OutputWriter     parent_type;

public:
	OutputWriterFile(void);
	OutputWriterFile(std::string input_file);
	OutputWriterFile(const self_type& other);
	~OutputWriterFile(void);

	bool open(const std::string& output_file);
	int writeHeaders(twk_header_type& twk_header);
	void writeFinal(void);
	void flush(void);
	void operator<<(buffer_type& buffer);
	void operator<<(const container_type& container){ parent_type::operator<<(container); }
	void writePrecompressedBlock(buffer_type& buffer, const U64& uncompressed_size);

private:
	void CheckOutputNames(const std::string& input);

public:
	std::string      filename;
	std::string      basePath;
	std::string      baseName;
	std::ofstream*   stream;
};

}
}

#endif /* IO_OUTPUT_WRITER_H_ */
