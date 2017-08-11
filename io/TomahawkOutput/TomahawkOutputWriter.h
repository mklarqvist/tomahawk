#ifndef TOMAHAWKOUTPUTWRITER_H_
#define TOMAHAWKOUTPUTWRITER_H_

#include "../BasicWriters.h"

namespace Tomahawk{
namespace IO {

class TomahawkOutputWriterInterface {
	typedef TomahawkOutputWriterInterface self_type;
	typedef GenericWriterInterace stream_type;
	typedef IO::WriterStandardOut cout_type;
	typedef IO::WriterFile file_type;

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
	virtual void operator<<(const void* entry) =0;
	virtual void write(const void* entry, const void* support) =0;
	inline void write(const char* data, const U32 length){ this->stream->write(&data[0], length); }

protected:
	std::string outFile;
	stream_type* stream;
};

// case binary
class TomahawkOutputWriter : public TomahawkOutputWriterInterface {
	typedef TomahawkOutputWriter self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::GZController tgzf_controller_type;
	typedef TomahawkOutputEntry entry_type;

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

	void write(const void* entry, const void* support){}

	void operator<<(const void* entry){
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
	typedef IO::BasicBuffer buffer_type;
	typedef IO::GZController tgzf_controller_type;
	typedef TomahawkOutputEntry entry_type;

public:
	TomahawkOutputWriterNatural(){}
	~TomahawkOutputWriterNatural(){
		// Flush upon termination
		this->flush();
		this->close();
	}

	void write(const void* entry, const void* support){
		*reinterpret_cast<std::ofstream*>(&this->stream->getStream()) << *reinterpret_cast<const entry_type*>(entry);
		this->stream->getStream() << '\n';
	}

	void operator<<(const void* entry){
		*reinterpret_cast<std::ofstream*>(&this->stream->getStream()) << *reinterpret_cast<const entry_type*>(entry);
		this->stream->getStream() << '\n';
	}
};


}
}

#endif /* TOMAHAWKOUTPUTWRITER_H_ */
