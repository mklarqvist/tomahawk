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
	virtual void write(const entry_type* const entry, const void* support) =0;
	inline void write(const char* data, const U32 length){ this->stream->write(&data[0], length); }

	template<class T> void write(const T& data){
		*reinterpret_cast<std::ofstream*>(&this->stream->getStream()) << data;
	}

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
		this->close();
		this->buffer.deleteAll();
	}

	inline void flush(void){ this->stream->flush(); }
	inline bool close(void){ this->stream->close(); return true; }

	void write(const entry_type* const entry, const void* support){}

	void operator<<(const entry_type* const entry){
		this->buffer.Add(reinterpret_cast<const char*>(entry), sizeof(entry_type));
		if(this->buffer.size() > this->flush_limit){
			this->controller.Deflate(this->buffer);
			this->stream->write(&this->controller.buffer_[0], this->controller.buffer_.size());
			this->controller.Clear();
			this->buffer.reset();
		}
	}

private:
	U32 flush_limit;
	buffer_type buffer;
	tgzf_controller_type controller;
};

// case binary
class TomahawkOutputWriterNatural : public TomahawkOutputWriterInterface {
	typedef TomahawkOutputWriterNatural self_type;

public:
	TomahawkOutputWriterNatural(){}
	~TomahawkOutputWriterNatural(){
		// Flush upon termination
		this->flush();
		this->close();
	}

	void write(const entry_type* const entry, const void* support){
		*reinterpret_cast<std::ofstream*>(&this->stream->getStream()) << *entry;
		this->stream->getStream() << '\n';
	}

	void operator<<(const entry_type* const entry){
		*reinterpret_cast<std::ofstream*>(&this->stream->getStream()) << *entry;
		this->stream->getStream() << '\n';
	}
};

}
}

#endif /* TOMAHAWKOUTPUTWRITER_H_ */
