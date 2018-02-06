#ifndef TOMAHAWK_TOMAHAWKOUTPUTMANAGER_H_
#define TOMAHAWK_TOMAHAWKOUTPUTMANAGER_H_

#include "../../io/BasicWriters.h"
#include "../../io/compression/TGZFController.h"
#include "../../support/MagicConstants.h"
#include "../../totempole/TotempoleContig.h"
#include "../../totempole/TotempoleMagic.h"
#include "../../totempole/TotempoleOutputEntry.h"
#include "TomahawkOutputEntry.h"
#include "TomahawkOutputLD.h"
#include "../../totempole/TotempoleReader.h"
#include "../../tomahawk/base/meta_entry.h"

#define SLAVE_FLUSH_LIMIT 10000000	// 10 MB

namespace Tomahawk{
namespace IO {

template <class T>
struct TomahawkOutputManager{
	typedef TomahawkOutputManager                     self_type;
	typedef IO::WriterFile                            writer_type;
	//typedef TomahawkBlock<const T, Tomahawk::Support::TomahawkRun<T>> controller_type;
	typedef Tomahawk::Support::TomahawkOutputLD       helper_type;
	typedef IO::BasicBuffer                           buffer_type;
	typedef TGZFController                            tgzf_controller;
	typedef IO::TomahawkOutputEntry                   entry_type;
	typedef Totempole::TotempoleOutputEntry           totempole_entry;
	typedef Totempole::TotempoleOutputEntryController totempole_controller_byte;
	typedef MetaEntry<T>                              meta_type;
	typedef Totempole::TotempoleEntry                 header_entry;

public:
	TomahawkOutputManager() :
		outCount(0),
		progressCount(0),
		totempole_blocks_written(0),
		writer(nullptr),
		writer_index(nullptr),
		buffer(2*SLAVE_FLUSH_LIMIT),
		sprintf_buffer(new char[255])
	{

	}

	~TomahawkOutputManager(){
		this->flushBlock();
		this->buffer.deleteAll();
		delete [] this->sprintf_buffer;
	}

	TomahawkOutputManager(const self_type& other) :
		outCount(0),
		progressCount(0),
		totempole_blocks_written(0),
		writer(other.writer),
		writer_index(other.writer_index),
		buffer(2*SLAVE_FLUSH_LIMIT),
		sprintf_buffer(new char[255])
	{
	}

	self_type& operator+=(const self_type& other){
		this->outCount += other.outCount;
		this->totempole_blocks_written += other.totempole_blocks_written;
		return(*this);
	}

	self_type& operator=(const U32 totempole_blocks){
		this->totempole_blocks_written = totempole_blocks;
		return(*this);
	}

	inline const U64& GetCounts(void) const{ return this->outCount; }
	inline void ResetProgress(void){ this->progressCount = 0; }
	inline const U32& GetProgressCounts(void) const{ return this->progressCount; }
	inline const U32& getTotempoleBlocks(void) const{ return(this->totempole_blocks_written); }

	bool Open(const std::string output, Totempole::TotempoleReader& totempole){
		if(output.size() == 0)
			return false;

		this->writer = new writer_type;
		this->writer_index = new writer_type;

		this->CheckOutputNames(output);
		this->filename = output;
		if(!this->writer->open(this->basePath + this->baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX)){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to open..." << std::endl;
			return false;
		}

		if(!this->writer_index->open(this->basePath + this->baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX)){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed open index..." << std::endl;
			return false;
		}

		if(!this->WriteHeader(totempole)){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to write header" << std::endl;
			return false;
		}

		return true;
	}

	void flushBlock(void){
		if(this->buffer.size() > 0){
			if(!this->compressor.Deflate(this->buffer)){
				std::cerr << Helpers::timestamp("ERROR","TGZF") << "Failed deflate DATA..." << std::endl;
				exit(1);
			}

			this->writer->getLock()->lock();
			this->entry.byte_offset = (U64)this->writer->getNativeStream().tellp();
			this->entry.uncompressed_size = this->buffer.size();
			this->writer->writeNoLock(this->compressor.buffer);
			this->entry.byte_offset_end = (U64)this->writer->getNativeStream().tellp();
			this->writer_index->getNativeStream() << this->entry;
			++this->totempole_blocks_written;
			//std::cerr << this->entry << std::endl;
			this->writer->getLock()->unlock();

			this->buffer.reset();
			this->compressor.Clear();
			this->entry.reset();
		}
	}

	bool finalise(void){
		// Make sure data is flushed
		this->writer->flush();
		this->writer_index->flush();

		// Update blocks written
		std::fstream re(this->basePath + this->baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX, std::ios::in | std::ios::out | std::ios::binary);
		if(!re.good()){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to reopen index..." << std::endl;
			return false;
		}

		re.seekg(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH + sizeof(float) + sizeof(U64) + sizeof(U32));
		if(!re.good()){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to seek in index..." << std::endl;
			return false;
		}

		re.write((char*)&this->totempole_blocks_written, sizeof(U32));
		if(!re.good()){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to update counts in index..." << std::endl;
			return false;
		}
		re.flush();
		re.close();

		return true;
	}

	void Add(const meta_type& meta_a, const meta_type& meta_b, const header_entry& header_a, const header_entry& header_b, const helper_type& helper){
		const U32 writePosA = meta_a.position << 2 | meta_a.phased << 1 | meta_a.missing;
		const U32 writePosB = meta_b.position << 2 | meta_b.phased << 1 | meta_b.missing;
		this->buffer += helper.controller;
		this->buffer += header_a.contigID;
		this->buffer += writePosA;
		this->buffer += header_b.contigID;
		this->buffer += writePosB;
		this->buffer << helper;
		++this->outCount;
		++this->progressCount;
		++this->entry.entries;

		if(this->buffer.size() > SLAVE_FLUSH_LIMIT)
			this->flushBlock();
	}

	void close(void){
		this->writer->flush();
		this->writer_index->flush();
		this->writer->close();
		this->writer_index->close();
		delete this->writer;
		delete this->writer_index;
	}

private:
	bool WriteHeader(Totempole::TotempoleReader& totempole){
		//typedef TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> header_type;
		std::ofstream& stream = this->writer->getNativeStream();
		std::ofstream& stream_index = this->writer_index->getNativeStream();

		TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> head(Tomahawk::Constants::WRITE_HEADER_LD_MAGIC, totempole.getSamples(), totempole.getContigs());
		TomahawkOutputSortHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> headIndex(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, totempole.getSamples(), totempole.getContigs());
		stream << head;
		stream_index << headIndex;

		// Write contig data to TWO
		// length | n_char | chars[0 .. n_char - 1]
		for(U32 i = 0; i < totempole.getContigs(); ++i)
			stream << *totempole.getContigBase(i);

		if(!totempole.writeLiterals(stream)){
			std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to write literals..." << std::endl;
			return false;
		}

		return(stream.good());
	}

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
	std::string     filename;
	std::string     basePath;
	std::string     baseName;
	U64             outCount;      // lines written
	U32             progressCount; // lines added since last flush
	U32             totempole_blocks_written;
	totempole_entry entry; // track stuff
	writer_type*    writer;	// writer
	writer_type*    writer_index; // writer index
	buffer_type     buffer; // internal buffer
	tgzf_controller compressor; // compressor
	char*           sprintf_buffer; // special buffer used for sprintf writing scientific output in natural mode
};

}
}

#endif /* TOMAHAWK_TOMAHAWKOUTPUTMANAGER_H_ */
