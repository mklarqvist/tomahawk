#ifndef TOMAHAWK_BASICWRITERS_H_
#define TOMAHAWK_BASICWRITERS_H_

#include <fstream>
#include <iostream>
#include <string>

#include "../support/helpers.h"
#include "../algorithm/spinlock.h"
#include "../io/BasicBuffer.h"
#include "../support/MagicConstants.h"

namespace Tomahawk {
namespace IO{

class GenericWriterInterace {
protected:
	typedef IO::BasicBuffer buffer_type;
	typedef Algorithm::SpinLock lock_type;

public:
	enum type {cout, file};
	enum compression {natural, binary};

public:
	GenericWriterInterace(){}
	virtual ~GenericWriterInterace(){}

	virtual bool open(void) =0;
	virtual bool open(const std::string output) =0;

	// Always the same but contents in buffer may be different
	virtual void operator<<(const buffer_type& buffer) =0;
	virtual void operator<<(void* entry) =0;

	virtual const size_t write(const char* data, const U64& length) =0;
	virtual std::ostream& getStream(void) =0;
	virtual void flush(void) =0;
	virtual bool close(void) =0;
	inline lock_type* getLock(void){ return(&this->lock); }

	virtual inline const size_t writeNoLock(const char* data, const U32 length) =0;
	virtual inline const size_t writeNoLock(const buffer_type& buffer) =0;

protected:
	lock_type lock;
};

class WriterStandardOut : public GenericWriterInterace{
public:
	WriterStandardOut(){}
	~WriterStandardOut(){
		// Flush upon termination
		this->flush();
	}

	bool open(void){ return true; }
	bool open(const std::string output){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Cannot set filename when destination is standard out..." << std::endl;
		return false;
	}
	void flush(void){ std::cout.flush(); }
	inline bool close(void){ return true; }
	inline std::ostream& getStream(void){ return(std::cout); }

	const size_t write(const char* data, const U64& length){
		this->lock.lock();
		std::cout.write(&data[0], length);
		this->lock.unlock();
		return(length);
	}

	inline const size_t writeNoLock(const char* data, const U32 length){
		std::cout.write(&data[0], length);
		return(length);
	}

	inline const size_t writeNoLock(const buffer_type& buffer){
		std::cout.write(&buffer.data[0], buffer.pointer);
		return(buffer.pointer);
	}

	void operator<<(void* entry){}

	void operator<<(const buffer_type& buffer){
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

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "WRITER") << "Opening output file: " << output << "..." << std::endl;

		return true;
	}

	inline std::ostream& getStream(void){ return(this->stream); }
	inline std::ofstream& getNativeStream(void){ return(this->stream); }
	inline void flush(void){ this->stream.flush(); }
	inline bool close(void){ this->stream.close(); return true; }

	void operator<<(const buffer_type& buffer){
		// Mutex lock; write; unlock
		// Note that this threads enter here at random
		// Extremely unlikely there is every any contention
		this->lock.lock();
		this->stream.write(&buffer.data[0], buffer.size());
		this->lock.unlock();
	}

	void operator<<(void* entry){}
	const size_t write(const char* data, const U64& length){
		this->lock.lock();
		this->stream.write(&data[0], length);
		this->lock.unlock();
		return(length);
	}

	inline const size_t writeNoLock(const char* data, const U32 length){
		this->stream.write(&data[0], length);
		return(length);
	}

	inline const size_t writeNoLock(const buffer_type& buffer){
		this->stream.write(&buffer.data[0], buffer.pointer);
		return(buffer.pointer);
	}

private:
	std::string outFile;
	std::ofstream stream;
};

} /* namespace IO */
} /* namespace Tomahawk */

#endif /* TOMAHAWK_BASICWRITERS_H_ */
