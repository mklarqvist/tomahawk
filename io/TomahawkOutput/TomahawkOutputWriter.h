#ifndef TOMAHAWKOUTPUTWRITER_H_
#define TOMAHAWKOUTPUTWRITER_H_

#include "../BasicWriters.h"

namespace Tomahawk{
namespace IO {

// case file
class TomahawkOutputWriter{
	typedef TomahawkOutputWriter self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::GZController tgzf_controller_type;
	typedef TomahawkOutputEntry entry_type;

public:
	TomahawkOutputWriter() : flush_limit(524288), stream(nullptr){}
	TomahawkOutputWriter(const U32 flush_limit) :
		flush_limit(flush_limit),
		stream(nullptr),
		buffer(flush_limit + 524288)
	{
		this->controller.buffer_.resize(this->buffer);
		delete this->stream;
	}

	~TomahawkOutputWriter(){
		// Flush upon termination
		this->flush();
		this->close();
		this->buffer.deleteAll();
	}

	bool open(void){
		this->stream = new IO::WriterStandardOut();
		if(this->stream->open()){
			return false;
		}
		return true;
	}

	bool open(const std::string output){
		if(output.size() == 0)
			return(this->open());

		this->stream = new IO::WriterFile();
		if(this->stream->open(output)){
			return false;
		}
		return true;
	}

	inline void flush(void){ this->stream->flush(); }
	inline bool close(void){ this->stream->close(); return true; }

	void operator<<(void* entry){
		this->buffer.Add(reinterpret_cast<const char*>(entry), sizeof(entry_type));
		if(this->buffer.size() > this->flush_limit){
			this->controller.Deflate(this->buffer);
			this->stream->write(&this->controller.buffer_[0], this->controller.buffer_.size());
			this->controller.Clear();
			this->buffer.reset();
		}
	}

	void write(const char* data, const U32 length){
		this->stream->write(&data[0], length);
	}

private:
	U32 flush_limit;
	std::string outFile;
	GenericWriterInterace* stream;
	buffer_type buffer;
	tgzf_controller_type controller;
};

}
}

#endif /* TOMAHAWKOUTPUTWRITER_H_ */
