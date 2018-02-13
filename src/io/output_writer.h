#ifndef IO_OUTPUT_WRITER_H_
#define IO_OUTPUT_WRITER_H_

#include "../support/MagicConstants.h"
#include "../support/simd_definitions.h"
#include "../support/helpers.h"
#include "../io/compression/TGZFController.h"
#include "../algorithm/spinlock.h"
#include "../index/tomahawk_header.h"
#include "../index/index_output_entry.h"
#include "../tomahawk/two/output_entry.h"
#include "../tomahawk/two/output_entry_support.h"
#include "../tomahawk/meta_entry.h"
#include "../index/index_entry.h"

namespace Tomahawk{
namespace IO{

class OutputWriter{
private:
	typedef OutputWriter                self_type;
	typedef TGZFController              compression_type;
	typedef Algorithm::SpinLock         spin_lock_type;
	typedef BasicBuffer                 buffer_type;
	typedef TomahawkHeader              twk_header_type;
	typedef Totempole::IndexOutputEntry index_entry_type;
	typedef OutputEntry                 entry_type;
	typedef Support::OutputEntrySupport entry_support_type;
	typedef Totempole::IndexEntry       header_entry_type;

public:
	OutputWriter(void) :
		owns_pointers(true),
		n_entries(0),
		n_progress_count(0),
		n_blocks(0),
		l_flush_limit(2000000),
		stream(nullptr),
		buffer(this->l_flush_limit*2),
		spin_lock(new spin_lock_type)
	{

	}

	OutputWriter(std::string input_file) :
		owns_pointers(true),
		n_entries(0),
		n_progress_count(0),
		n_blocks(0),
		l_flush_limit(2000000),
		stream(new std::ofstream(input_file, std::ios::binary | std::ios::out)),
		buffer(this->l_flush_limit*2),
		spin_lock(new spin_lock_type)
	{

	}

	OutputWriter(const self_type& other) :
		owns_pointers(false),
		n_entries(other.n_entries),
		n_progress_count(other.n_progress_count),
		n_blocks(other.n_blocks),
		l_flush_limit(other.l_flush_limit),
		stream(other.stream),
		buffer(other.buffer.capacity()),
		spin_lock(other.spin_lock)
	{

	}

	~OutputWriter(void){
		if(this->owns_pointers){
			this->stream->flush();
			this->stream->close();
			delete this->stream;
			delete this->spin_lock;
		}
	}

	bool open(const std::string& output_file){
		if(output_file.size() == 0)
			return false;

		this->CheckOutputNames(output_file);
		this->filename = output_file;

		this->stream = new std::ofstream(this->basePath + this->baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX, std::ios::binary | std::ios::out);
		if(this->stream->good() == false){
			std::cerr << "Failed to open: " << output_file << std::endl;
			return false;
		}

		return true;
	}

	int WriteHeaders(twk_header_type& twk_header){
		const std::string command = "##tomahawk_calcCommand=" + Helpers::program_string();
		twk_header.getLiterals() += command;
		// Set file type to TWO
		twk_header.magic_.file_type = 1;

		return(twk_header.write(*this->stream));
	}

	void flush(void){
		if(this->buffer.size() > 0){
			if(!this->compressor.Deflate(this->buffer)){
				std::cerr << Helpers::timestamp("ERROR","TGZF") << "Failed deflate DATA..." << std::endl;
				exit(1);
			}

			this->spin_lock->lock();
			this->index_entry.byte_offset_from = (U64)this->stream->tellp();
			this->index_entry.uncompressed_size = this->buffer.size();
			this->stream->write(this->compressor.buffer.data(), this->compressor.buffer.size());
			this->index_entry.byte_offset_to = (U64)this->stream->tellp();
			//*this->stream << this->index_entry;
			//std::cerr << this->index_entry.byte_offset_from << "->" << this->index_entry.byte_offset_to << " for " << this->index_entry.n_entries << " of " << this->index_entry.uncompressed_size << std::endl;
			++this->n_blocks;

			this->spin_lock->unlock();

			this->buffer.reset();
			this->compressor.Clear();
			this->index_entry.reset();
		}
	}

	inline self_type& operator+=(const self_type& other){
		this->n_entries += other.n_entries;
		this->n_blocks  += other.n_blocks;
		return(*this);
	}

	inline self_type& operator=(const self_type& other){
		this->n_blocks  = other.n_blocks;
		this->n_entries = other.n_entries;
		return(*this);
	}

	template <class T>
	void Add(const MetaEntry<T>& meta_a, const MetaEntry<T>& meta_b, const header_entry_type& header_a, const header_entry_type& header_b, const entry_support_type& helper){
		const U32 writePosA = meta_a.position << 2 | meta_a.phased << 1 | meta_a.missing;
		const U32 writePosB = meta_b.position << 2 | meta_b.phased << 1 | meta_b.missing;
		this->buffer += helper.controller;
		this->buffer += header_a.contigID;
		this->buffer += writePosA;
		this->buffer += header_b.contigID;
		this->buffer += writePosB;
		this->buffer << helper;
		++this->n_entries;
		++this->n_progress_count;
		++this->index_entry.n_entries;

		if(this->buffer.size() > this->l_flush_limit)
			this->flush();
	}

	inline void ResetProgress(void){ this->n_progress_count = 0; }
	inline const U32& getProgressCounts(void) const{ return this->n_progress_count; }

	inline void operator<<(const entry_type& entry){
		this->buffer << entry;
		++this->n_entries;
	}

	void operator<<(buffer_type& buffer){
		if(!this->compressor.Deflate(buffer)){
			std::cerr << "failed to add" << std::endl;
			return;
		}
		this->spin_lock->lock();
		this->stream->write(this->compressor.buffer.data(), this->compressor.buffer.size());
		this->stream->flush();
		this->spin_lock->unlock();
		this->compressor.buffer.reset();
		buffer.reset();
	}

private:
	void CheckOutputNames(const std::string& input){
		std::vector<std::string> paths = Helpers::filePathBaseExtension(input);
		this->basePath = paths[0];
		if(this->basePath.size() > 0)
			this->basePath += '/';

		if(paths[3].size() == Tomahawk::Constants::OUTPUT_LD_SUFFIX.size() && strncasecmp(&paths[3][0], &Tomahawk::Constants::OUTPUT_LD_SUFFIX[0], Tomahawk::Constants::OUTPUT_LD_SUFFIX.size()) == 0)
			this->baseName = paths[2];
		else this->baseName = paths[1];
	}

private:
	std::string       filename;
	std::string       basePath;
	std::string       baseName;
	bool              owns_pointers;
	U64               n_entries;        // number of entries written
	U32               n_progress_count; // lines added since last flush
	U32               n_blocks;         // number of index blocks writtenflush_limit
	U32               l_flush_limit;
	index_entry_type  index_entry;      // keep track of sort order
	std::ofstream*    stream;
	buffer_type       buffer;
	compression_type  compressor;
	spin_lock_type*   spin_lock;
};

}
}

#endif /* IO_OUTPUT_WRITER_H_ */