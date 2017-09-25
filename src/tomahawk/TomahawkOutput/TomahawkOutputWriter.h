#ifndef TOMAHAWKOUTPUTWRITER_H_
#define TOMAHAWKOUTPUTWRITER_H_

#include "../../io/BasicWriters.h"
#include "../../totempole/TotempoleMagic.h"
#include "../../totempole/TotempoleOutputEntry.h"

namespace Tomahawk{
namespace IO {

class TomahawkOutputWriterInterface {
private:
	typedef TomahawkOutputWriterInterface self_type;

protected:
	typedef GenericWriterInterace stream_type;
	typedef IO::WriterStandardOut cout_type;
	typedef IO::WriterFile file_type;
	typedef IO::BasicBuffer buffer_type;

protected:
	typedef TomahawkOutputEntry entry_type;
	typedef TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> header_type;
	typedef Totempole::TotempoleContigBase contig_type;

public:
	typedef Algorithm::SpinLock lock_type;

public:
	TomahawkOutputWriterInterface(const contig_type* contigs, const header_type* header) :
		stream(nullptr),
		header(header),
		contigs(contigs)
	{}
	virtual ~TomahawkOutputWriterInterface(){ delete this->stream; }

	virtual bool open(void){
		this->stream = new cout_type();
		if(!this->stream->open())
			return false;

		return true;
	}

	virtual bool open(const std::string output){
		if(output.size() == 0)
			return(this->open());

		this->outFile = output;
		this->CheckOutputNames(output);

		this->stream = new file_type();
		if(!this->stream->open(this->basePath + this->baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX)){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to open..." << std::endl;
			return false;
		}

		return true;
	}

	virtual inline void flush(void){ this->stream->flush(); }
	virtual inline bool close(void){ this->stream->close(); return true; }
	virtual void operator<<(const entry_type* const entryentry) =0;
	virtual void operator<<(const entry_type& entryentry) =0;
	virtual inline void write(const char* data, const U32 length){ this->stream->write(&data[0], length); }
	virtual inline void writeNoLock(const char* data, const U32 length){ this->stream->writeNoLock(&data[0], length); }
	virtual inline const size_t write(buffer_type& buffer){ return(this->stream->write(&buffer.data[0], buffer.size())); }
	virtual inline const size_t writeNoLock(buffer_type& buffer){ return(this->stream->writeNoLock(&buffer.data[0], buffer.size())); }
	inline stream_type* getStream(void){ return(this->stream); }
	inline lock_type* getLock(void){ return(this->stream->getLock()); }

	template<class T>
	void write(const T& data){
		*reinterpret_cast<std::ofstream*>(&this->stream->getStream()) << data;
	}

	virtual void writeHeader(std::string& literals) =0;
	virtual void writeEOF(void) =0;

	void CheckOutputNames(const std::string& input){
		std::vector<std::string> paths = Helpers::filePathBaseExtension(input);
		this->basePath = paths[0];
		if(this->basePath.size() > 0)
			this->basePath += '/';

		if(paths[3].size() == Tomahawk::Constants::OUTPUT_LD_SUFFIX.size() && strncasecmp(&paths[3][0], &Tomahawk::Constants::OUTPUT_LD_SUFFIX[0], Tomahawk::Constants::OUTPUT_LD_SUFFIX.size()) == 0)
			this->baseName = paths[2];
		else this->baseName = paths[1];
	}

protected:
	std::string outFile;
	std::string basePath;
	std::string baseName;

	stream_type* stream;
	const header_type* header;
	const contig_type* contigs;
};

// case binary
class TomahawkOutputWriter : public TomahawkOutputWriterInterface {
private:
	typedef TomahawkOutputWriter self_type;

protected:
	typedef IO::BasicBuffer buffer_type;
	typedef IO::TGZFController tgzf_controller_type;

public:
	TomahawkOutputWriter(const contig_type* contigs, const header_type* header) : TomahawkOutputWriterInterface(contigs, header), flush_limit(524288){}
	TomahawkOutputWriter(const contig_type* contigs, const header_type* header, const U32 flush_limit) :
		TomahawkOutputWriterInterface(contigs, header),
		flush_limit(flush_limit),
		buffer(flush_limit + 2048)
	{
		this->controller.buffer.resize(this->buffer);
	}

	// if binary
	// does not support cout writer
	// also open two.twi
	// if not ending in .two add .two

	virtual ~TomahawkOutputWriter(){
		// Flush upon termination
		this->flush();
		this->writeEOF();
		this->close();
		this->buffer.deleteAll();
	}

	inline void flush(void){
		if(this->buffer.size() > 0){
			this->controller.Deflate(this->buffer);
			this->stream->write(&this->controller.buffer[0], this->controller.buffer.size());
			this->controller.Clear();
			this->buffer.reset();
		}
		this->stream->flush();
	}
	virtual inline bool close(void){ this->stream->close(); return true; }

	virtual void operator<<(const entry_type* const entry){
		this->buffer.Add(reinterpret_cast<const char*>(entry), sizeof(entry_type));
		if(this->buffer.size() > this->flush_limit){
			this->controller.Deflate(this->buffer);
			this->stream->write(&this->controller.buffer[0], this->controller.buffer.size());
			this->controller.Clear();
			this->buffer.reset();
		}
	}

	virtual void operator<<(const entry_type& entry){
		this->buffer.Add(reinterpret_cast<const char*>(&entry), sizeof(entry_type));
		if(this->buffer.size() > this->flush_limit){
			this->controller.Deflate(this->buffer);
			this->stream->write(&this->controller.buffer[0], this->controller.buffer.size());
			this->controller.Clear();
			this->buffer.reset();
		}
	}

	virtual inline const size_t write(buffer_type& buffer){
		this->controller.Clear();
		this->controller.Deflate(buffer);
		this->stream->write(&this->controller.buffer[0], this->controller.buffer.size());
		return(this->controller.buffer.size());
	}

	virtual inline const size_t writeNoLock(buffer_type& buffer){
		this->controller.Clear();
		this->controller.Deflate(buffer);
		this->stream->writeNoLock(&this->controller.buffer[0], this->controller.buffer.size());
		return(this->controller.buffer.size());
	}

	void writeHeader(std::string& literals){
		std::ofstream& __stream = *reinterpret_cast<std::ofstream*>(&this->stream->getStream());
		__stream << *this->header;
		for(U32 i = 0; i < this->header->n_contig; ++i)
			__stream << this->contigs[i];

		buffer_type bufferInternal(&literals[0], literals.size());
		if(!this->controller.Deflate(bufferInternal)){
			std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to deflate!" << std::endl;
			return;
		}

		__stream.write(&this->controller.buffer.data[0], this->controller.buffer.pointer);
		this->controller.Clear();
	};

	U64 getUncompressedSize(void) const{ return(this->controller.buffer.size()); }

	// There is no EOF
	void writeEOF(void){};

protected:
	U32 flush_limit;
	buffer_type buffer;
	tgzf_controller_type controller;
};

// case binary + index
class TomahawkOutputWriterIndex : public TomahawkOutputWriter{
	typedef TomahawkOutputWriter parent_type;
	typedef Totempole::TotempoleOutputEntry totempole_type;
	typedef Tomahawk::IO::TomahawkOutputSortHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> toi_header_type;
	typedef Totempole::TotempoleOutputSortedIndex index_type;

public:
	TomahawkOutputWriterIndex(const contig_type* contigs, const header_type* header, const toi_header_type& toi_header) :
		TomahawkOutputWriter(contigs, header),
		current_blockID(0),
		toi_header(toi_header),
		index(header->n_contig, contigs)
	{

	}

	TomahawkOutputWriterIndex(const contig_type* contigs, const header_type* header) :
		TomahawkOutputWriter(contigs, header),
		current_blockID(0),
		index(header->n_contig, contigs)
	{
		this->toi_header = toi_header_type(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, header->samples, header->n_contig);
	}

	TomahawkOutputWriterIndex(const contig_type* contigs, const header_type* header, const U32 flush_limit) :
		TomahawkOutputWriter(contigs, header, flush_limit),
		current_blockID(0),
		index(header->n_contig, contigs)
	{
		this->toi_header = toi_header_type(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, header->samples, header->n_contig);
	}

	~TomahawkOutputWriterIndex(){}

	bool open(const std::string output){
		if(output.size() == 0){
			std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Writing to cout is illegal..." << std::endl;
			return false;
		}

		this->outFile = output;
		this->CheckOutputNames(output);

		this->stream = new file_type();
		if(!this->stream->open(this->basePath + this->baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX)){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to open..." << std::endl;
			return false;
		}

		if(!this->stream_index.open(this->basePath + this->baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX)){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed open index..." << std::endl;
			return false;
		}

		return true;
	}

	inline bool close(void){
		this->stream->flush();
		this->stream->close();
		this->stream_index.flush();
		this->stream_index.close();
		return true;
	}

	inline void flush(void){
		this->flushBlock();
		this->stream_index.flush();
		this->stream->flush();
	}

	void operator<<(const entry_type* const entry){
		++this->totempole_entry.entries;

		this->index.update(*entry, this->current_blockID, this->buffer.size() / sizeof(entry_type));
		this->buffer.Add(reinterpret_cast<const char* const>(&entry), sizeof(entry_type));

		if(this->buffer.size() > this->flush_limit)
			this->flushBlock();
	}

	void operator<<(const entry_type& entry){
		++this->totempole_entry.entries;

		this->index.update(entry, this->current_blockID, this->buffer.size() / sizeof(entry_type));
		this->buffer.Add(reinterpret_cast<const char* const>(&entry), sizeof(entry_type));

		if(this->buffer.size() > this->flush_limit)
			this->flushBlock();
	}

	void writeHeader(std::string& literals){
		parent_type::writeHeader(literals);
		this->stream_index.getNativeStream() << toi_header;
	};

	void flushBlock(void){
		if(this->buffer.size() == 0)
			return;

		IO::WriterFile& stream = *reinterpret_cast<IO::WriterFile*>(this->stream);
		this->controller.Deflate(this->buffer);
		this->totempole_entry.byte_offset = stream.getNativeStream().tellp();
		this->stream->write(&this->controller.buffer[0], this->controller.buffer.size());
		this->totempole_entry.byte_offset_end = stream.getNativeStream().tellp();
		this->totempole_entry.uncompressed_size = this->controller.buffer.size();
		this->controller.Clear();
		this->buffer.reset();
		this->stream_index.getNativeStream() << this->totempole_entry;
		this->totempole_entry.reset();
		++this->current_blockID;
	}

	index_type& getIndex(void){ return(this->index); }

	bool finalize(bool output_index){
		// Totempole Output Index
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

		re.write((char*)&current_blockID, sizeof(U32));
		if(!re.good()){
			std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to update counts in index..." << std::endl;
			return false;
		}
		re.close();

		// Write sorted index
		if(output_index)
			this->stream_index.getNativeStream() << this->index;

		return true;
	}

private:
	U32 current_blockID;
	toi_header_type toi_header;
	totempole_type totempole_entry;
	entry_type entry_prev;
	file_type stream_index;
	index_type index;
};

// case natural
class TomahawkOutputWriterNatural : public TomahawkOutputWriterInterface {
	typedef TomahawkOutputWriterNatural self_type;

public:
	TomahawkOutputWriterNatural(const contig_type* contigs, const header_type* header) : TomahawkOutputWriterInterface(contigs, header){}
	~TomahawkOutputWriterNatural(){
		// Flush upon termination
		this->flush();
		this->writeEOF();
		this->close();
	}

	void operator<<(const entry_type* const entry){
		entry->write(*reinterpret_cast<std::ofstream*>(&this->stream->getStream()), this->contigs);
	}

	void operator<<(const entry_type& entry){
		entry.write(*reinterpret_cast<std::ofstream*>(&this->stream->getStream()), this->contigs);
	}

	void writeHeader(std::string& literals){
		const std::string header = "FLAG\tAcontigID\tAposition\tBcontigID\tBpositionID\tp1\tp2\tq1\tq2\tD\tDprime\tR2\tP\tchiSqFisher\tchiSqModel\n";
		if(literals.size() > 0)
			this->stream->getStream() << literals << '\n' << header;
		else
			this->stream->getStream() << header;
	};

	// There is no EOF
	void writeEOF(void){};
};

}
}

#endif /* TOMAHAWKOUTPUTWRITER_H_ */
