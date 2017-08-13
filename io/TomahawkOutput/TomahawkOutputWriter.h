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
	typedef Totempole::TotempoleContigBase contig_type;

public:
	TomahawkOutputWriterInterface(const contig_type* contigs, const header_type* header) :
		stream(nullptr),
		header(header),
		contigs(contigs)
	{}
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
	const header_type* header;
	const contig_type* contigs;
};

// case binary
class TomahawkOutputWriter : public TomahawkOutputWriterInterface {
	typedef TomahawkOutputWriter self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::TGZFController tgzf_controller_type;

public:
	TomahawkOutputWriter(const contig_type* contigs, const header_type* header) : TomahawkOutputWriterInterface(contigs, header), flush_limit(524288){}
	TomahawkOutputWriter(const contig_type* contigs, const header_type* header, const U32 flush_limit) :
		TomahawkOutputWriterInterface(contigs, header),
		flush_limit(flush_limit),
		buffer(flush_limit + 524288)
	{
		this->controller.buffer.resize(this->buffer);
	}

	~TomahawkOutputWriter(){
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
	inline bool close(void){ this->stream->close(); return true; }

	void operator<<(const entry_type* const entry){
		this->buffer.Add(reinterpret_cast<const char*>(entry), sizeof(entry_type));
		if(this->buffer.size() > this->flush_limit){
			this->controller.Deflate(this->buffer);
			this->stream->write(&this->controller.buffer[0], this->controller.buffer.size());
			this->controller.Clear();
			this->buffer.reset();
		}
	}

	void writeHeader(void){
		*reinterpret_cast<std::ofstream*>(&this->stream->getStream()) << *this->header;
		for(U32 i = 0; i < this->header->n_contig; ++i)
			*reinterpret_cast<std::ofstream*>(&this->stream->getStream()) << this->contigs[i];

	};

	// There is no EOF
	void writeEOF(void){};

private:
	U32 flush_limit;
	buffer_type buffer;
	tgzf_controller_type controller;
};

// case binary
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

	void writeHeader(void){
		const std::string header = "FLAG\tAcontigID\tAposition\tBcontigID\tBpositionID\tp1\tp2\tq1\tq2\tD\tDprime\tR2\tP\tchiSqFisher\tchiSqModel\n";
		this->stream->getStream() << header;
	};

	// There is no EOF
	void writeEOF(void){};

};

}
}

#endif /* TOMAHAWKOUTPUTWRITER_H_ */
