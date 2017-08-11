#ifndef TOMAHAWKOUTPUTWRITER_H_
#define TOMAHAWKOUTPUTWRITER_H_

namespace Tomahawk{
namespace IO {

// case file
class TomahawkOutputWriter : public GenericWriterInterace{
	typedef TomahawkOutputWriter self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::GZController tgzf_controller_type;
	typedef TomahawkOutputEntry entry_type;

public:
	TomahawkOutputWriter() : flush_limit(0){}
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

	bool open(void){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "No output name provided..." << std::endl;
		return false;
	}

	bool open(const std::string output){
		if(output.length() == 0){
			std::cerr << Helpers::timestamp("ERROR", "WRITER") << "No output name provided..." << std::endl;
			return false;
		}

		this->stream.open(output, std::ios::binary | std::ios::out);
		if(!this->stream.good()){
			std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not open output file: " << output << "..." << std::endl;
			return false;
		}

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "WRITER") << "Opening output file: " << output << "..." << std::endl;

		return true;
	}

	inline std::ostream& getStream(void){ return(this->stream); }
	inline std::ofstream& getNativeStream(void){ return(this->stream); }
	inline void flush(void){ this->stream.flush(); }
	inline bool close(void){ this->stream.close(); return true; }

	void operator<<(void* entry){
		const entry_type* const e = (const entry_type* const)entry;
		//std::cerr << "adding: " << *e << std::endl;
		this->buffer << *e;
		if(this->buffer.size() > this->flush_limit){
			std::cerr << this->buffer.size() << ">" << this->flush_limit << std::endl;
			this->buffer.reset();
		}
	}

	void operator<<(const buffer_type& buffer){
		//this->stream.write(&buffer.data[0], buffer.size());
	}

	void write(const char* data, const U32 length){
		this->stream.write(&data[0], length);
	}

private:
	U32 flush_limit;
	std::string outFile;
	std::ofstream stream;
	buffer_type buffer;
	tgzf_controller_type controller;
};

}
}

#endif /* TOMAHAWKOUTPUTWRITER_H_ */
