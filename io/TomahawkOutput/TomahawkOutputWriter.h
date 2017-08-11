#ifndef TOMAHAWKOUTPUTWRITER_H_
#define TOMAHAWKOUTPUTWRITER_H_

#include "../BasicWriters.h"
#include "../../totempole/TotempoleMagic.h"

namespace Tomahawk{
namespace IO {

class TomahawkOutputWriterInterface {
	typedef TomahawkOutputWriterInterface self_type;
	typedef GenericWriterInterace stream_type;
	typedef IO::WriterStandardOut cout_type;
	typedef IO::WriterFile file_type;

protected:
	typedef TomahawkOutputEntry entry_type;
	typedef TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> header_type;

public:
	TomahawkOutputWriterInterface() : stream(nullptr){}
	virtual ~TomahawkOutputWriterInterface(){ delete this->stream; }

	bool open(void){
		this->stream = new cout_type();
		if(this->stream->open()){
			return false;
		}

		return true;
	}

	bool open(const std::string output){
		if(output.size() == 0)
			return(this->open());

		this->outFile = output;
		this->stream = new file_type();
		if(this->stream->open(output)){
			return false;
		}

		return true;
	}

	virtual inline void flush(void){ this->stream->flush(); }
	virtual inline bool close(void){ this->stream->close(); return true; }
	virtual void operator<<(const entry_type* const entryentry) =0;
	inline void write(const char* data, const U32 length){ this->stream->write(&data[0], length); }

	template<class T> void write(const T& data){
		*reinterpret_cast<std::ofstream*>(&this->stream->getStream()) << data;
	}

	virtual void writeHeader(void) =0;
	virtual void writeEOF(void) =0;

protected:
	std::string outFile;
	stream_type* stream;
};

// case binary
class TomahawkOutputWriter : public TomahawkOutputWriterInterface {
	typedef TomahawkOutputWriter self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::GZController tgzf_controller_type;

public:
	TomahawkOutputWriter() : flush_limit(524288){}
	TomahawkOutputWriter(const U32 flush_limit) :
		flush_limit(flush_limit),
		buffer(flush_limit + 524288)
	{
		this->controller.buffer_.resize(this->buffer);
	}

	~TomahawkOutputWriter(){
		// Flush upon termination
		this->flush();
		this->writeEOF();
		this->close();
		this->buffer.deleteAll();
	}

	inline void flush(void){ this->stream->flush(); }
	inline bool close(void){ this->stream->close(); return true; }

	void operator<<(const entry_type* const entry){
		this->buffer.Add(reinterpret_cast<const char*>(entry), sizeof(entry_type));
		if(this->buffer.size() > this->flush_limit){
			this->controller.Deflate(this->buffer);
			this->stream->write(&this->controller.buffer_[0], this->controller.buffer_.size());
			this->controller.Clear();
			this->buffer.reset();
		}
	}

	void writeHeader(void){	};
	void writeEOF(void){};

private:
	U32 flush_limit;
	buffer_type buffer;
	tgzf_controller_type controller;
};

// case binary
class TomahawkOutputWriterNatural : public TomahawkOutputWriterInterface {
	typedef TomahawkOutputWriterNatural self_type;
	typedef Totempole::TotempoleContigBase contig_type;

public:
	TomahawkOutputWriterNatural() : contigs(nullptr){}
	~TomahawkOutputWriterNatural(){
		// Flush upon termination
		this->flush();
		this->writeEOF();
		this->close();
	}

	void operator<<(const entry_type* const entry){
		entry->write(*reinterpret_cast<std::ofstream*>(&this->stream->getStream()), this->contigs);
	}

	void setContigs(const contig_type* const contigs){ this->contigs = contigs; }
	void writeHeader(void){
		const std::string header = "FLAG\tAcontigID\tAposition\tBcontigID\tBpositionID\tp1\tp2\tq1\tq2\tD\tDprime\tR2\tP\tchiSqFisher\tchiSqModel\n";
		this->stream->getStream() << header;
	};
	void writeEOF(void){};

private:
	const contig_type* contigs;
};

}
}

#endif /* TOMAHAWKOUTPUTWRITER_H_ */
