#ifndef TOMAHAWK_BASICWRITERS_H_
#define TOMAHAWK_BASICWRITERS_H_

#include <fstream>
#include <iostream>
#include <string>

#include "../helpers.h"
#include "../algorithm/spinlock.h"
#include "../io/BasicBuffer.h"
#include "TomahawkOutputLD.h"
#include "../totempole/TotempoleReader.h"

namespace Tomahawk {
namespace IO{

// temp
class GenericWriterInterace{
protected:
	typedef Support::TomahawkOutputLD helper_type;
	typedef TotempoleReader totempole_type;

public:
	enum type {cout, file};
	enum compression {natural, binary};

public:
	GenericWriterInterace(){}
	virtual ~GenericWriterInterace(){}

	virtual bool open(void) =0;
	virtual bool open(const std::string output) =0;

	// Always the same but contents in buffer may be different
	virtual void operator<<(const IO::BasicBuffer& buffer) =0;

	// Header output
	virtual bool writeHeader(void) =0;
	virtual bool writeHeader(const totempole_type& totempole) =0;

	virtual void flush(void) =0;
	virtual bool close(void) =0;

protected:
	Algorithm::SpinLock lock;
};

class WriterStandardOut : public GenericWriterInterace{
public:
	WriterStandardOut(){}
	~WriterStandardOut(){
		// Flush upon termination
		this->flush();
		this->close();
	}

	bool open(void){ return true; }
	bool open(const std::string output){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Cannot set output filename when destination is standard out..." << std::endl;
		return false;
	}
	void flush(void){ std::cout.flush(); }
	inline bool close(void){ return true; }

	bool writeHeader(void){ return false; }
	bool writeHeader(const totempole_type& totempole){ return false; }

	inline void operator<<(const IO::BasicBuffer& buffer){
		// Mutex lock; write; unlock
		// Note that this threads enter here at random
		// Extremely unlikely there is every any contention
		this->lock.lock();
		std::cout.write(&buffer.data[0], buffer.size());
		this->lock.unlock();
	}
};

// case file
class WriterFile : public GenericWriterInterace{
public:
	WriterFile(){}
	~WriterFile(){
		// Flush upon termination
		this->flush();
		this->close();
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

		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Opening output file: " << output << "..." << std::endl;

		return true;
	}

	bool writeHeader(void){ return false; }
	bool writeHeader(const totempole_type& totempole){ return false; }

	inline void flush(void){ this->stream.flush(); }
	inline bool close(void){ this->stream.close(); return true; }

	inline void operator<<(const IO::BasicBuffer& buffer){
		// Mutex lock; write; unlock
		// Note that this threads enter here at random
		// Extremely unlikely there is every any contention
		this->lock.lock();
		this->stream.write(&buffer.data[0], buffer.size());
		this->lock.unlock();
	}

private:
	std::string outFile;
	std::ofstream stream;
};
////////////

} /* namespace IO */
} /* namespace Tomahawk */

#endif /* TOMAHAWK_BASICWRITERS_H_ */
