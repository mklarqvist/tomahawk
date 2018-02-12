#ifndef IO_OUTPUT_WRITER_H_
#define IO_OUTPUT_WRITER_H_

#include "../support/MagicConstants.h"
#include "../support/simd_definitions.h"
#include "../support/helpers.h"

namespace Tomahawk{
namespace IO{

class OutputWriter{
private:
	typedef OutputWriter        self_type;
	typedef TGZFController      compression_type;
	typedef Algorithm::SpinLock spin_lock_type;
	typedef BasicBuffer         buffer_type;
	typedef TomahawkHeader      twk_header_type;

public:
	OutputWriter(void) :
		owns_pointers(true),
		stream(nullptr),
		spin_lock(new spin_lock_type)
	{

	}

	OutputWriter(std::string input_file) :
		owns_pointers(true),
		stream(new std::ofstream(input_file, std::ios::binary | std::ios::out)),
		spin_lock(new spin_lock_type)
	{

	}

	OutputWriter(const self_type& other) :
		owns_pointers(false),
		stream(other.stream),
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

	bool open(const std::string& input_file){
		this->stream = new std::ofstream(input_file, std::ios::binary | std::ios::out);
		if(this->stream->good() == false){
			std::cerr << "faailed to open: " << input_file << std::endl;
			return false;
		}

		return true;
	}

	int WriteHeaders(twk_header_type& twk_header){
		const std::string command = "##tomahawk_importCommand=" + Helpers::program_string();;
		twk_header.getLiterals() += command;

		return(twk_header.write(*this->stream));
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
	}

private:
	bool              owns_pointers;
	std::ofstream*    stream;
	compression_type  compressor;
	spin_lock_type*   spin_lock;
};

}
}



#endif /* IO_OUTPUT_WRITER_H_ */
