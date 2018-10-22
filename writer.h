#ifndef TWK_WRITER_H_
#define TWK_WRITER_H_

#include <ostream>
#include "index.h"
#include "spinlock.h"
#include "core.h"

namespace tomahawk {

struct twk_writer_t {

	twk_writer_t() : owner(false), stream(nullptr){}
	virtual ~twk_writer_t(){ if(owner) delete stream; }

	virtual bool Open(const std::string& file) =0;

	void Add(twk_buffer_t& buffer){
		spinlock.lock();
		stream->write(buffer.data(), buffer.size());
		spinlock.unlock();
		buffer.reset();
	}

	void AddTwoBlock(twk_buffer_t& buffer){
		spinlock.lock();
		uint8_t marker = 1;
		stream->write(reinterpret_cast<const char*>(&marker), sizeof(uint8_t));
		stream->write(buffer.data(), buffer.size());
		spinlock.unlock();
		buffer.reset();
	}

	template <class T> twk_writer_t& operator<<(const T& val){ *this->stream << val; return(*this); }

	SpinLock spinlock;
	IndexEntry ientry;
	Index index;
	bool owner;
	std::ostream* stream;
};

struct twk_writer_stream : public twk_writer_t {
	inline bool Open(const std::string& file){
		stream = &std::cout;
		owner = false;
		return true;
	}
};

struct twk_writer_file : public twk_writer_t {
	inline bool Open(const std::string& file){
		owner = true;
		stream = new std::ofstream;
		std::ofstream* outstream = reinterpret_cast<std::ofstream*>(stream);
		outstream->open(file,std::ios::out | std::ios::binary);
		if(!outstream->good()){
			std::cerr << "failed to open" << std::endl;
			return false;
		}
		//
		return true;
	}
};


}

#endif /* TWK_WRITER_H_ */
