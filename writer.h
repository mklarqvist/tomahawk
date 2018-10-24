#ifndef TWK_WRITER_H_
#define TWK_WRITER_H_

#include <ostream>
#include "index.h"
#include "spinlock.h"
#include "core.h"

namespace tomahawk {

struct twk_writer_t {

	twk_writer_t() : buf(nullptr), stream(nullptr){}
	virtual ~twk_writer_t(){ stream.flush();  }

	virtual bool Open(const std::string& file) =0;

	void Add(twk_buffer_t& buffer){
		spinlock.lock();
		stream.write(buffer.data(), buffer.size());
		spinlock.unlock();
		buffer.reset();
	}

	/**<
	 * Write the equivalent of a twk_oblock_two_t block. Takes the arguments
	 * for the uncompressed byte size and the compressed byte size together
	 * with the actual byte buffer.
	 * @param b_unc  Uncompresed size in bytes.
	 * @param b_comp Compresed size in bytes.
	 * @param obuf   Src buffer.
	 */
	void Add(const uint32_t b_unc, const uint32_t b_comp, twk_buffer_t& obuf){
		spinlock.lock();
		uint8_t marker = 1; // non-termination marker
		SerializePrimitive(marker, stream);
		SerializePrimitive(b_unc, stream); // uncompressed size
		SerializePrimitive(b_comp, stream); // compressed size
		stream.write(obuf.data(), obuf.size()); // actual data
		spinlock.unlock();
		obuf.reset();
	}

	inline void write(const char* data, const uint32_t l){ stream.write(data, l); }
	inline void flush(){ stream.flush(); }
	inline bool good() const{ return(stream.good()); }
	inline std::streambuf* rdbuf(){ return(buf); }
	virtual void close() =0;
	virtual bool is_open() const =0;
	inline int64_t tellp(){ return(stream.tellp()); }
	inline void seekp(const uint64_t p){ stream.seekp(p); }

	template <class T> twk_writer_t& operator<<(const T& val){ this->stream << val; return(*this); }

public:
	std::streambuf* buf;
	SpinLock spinlock;
	IndexEntry ientry;
	Index index;
	std::ostream stream;
};

struct twk_writer_stream : public twk_writer_t {
	twk_writer_stream(){
		buf = std::cout.rdbuf();
		stream.basic_ios<char>::rdbuf(buf);
	}

	bool Open(const std::string& file){
		buf = std::cout.rdbuf();
		stream.basic_ios<char>::rdbuf(buf);
		return true;
	}

	void close(){}
	bool is_open() const{ return true; }

};

struct twk_writer_file : public twk_writer_t {
	bool Open(const std::string& file){
		out.open(file,std::ios::out | std::ios::binary);
		if(!out.good()){
			std::cerr << "failed to open" << std::endl;
			return false;
		}
		buf = out.rdbuf();
		stream.rdbuf(buf);
		return true;
	}

	void close(){ out.close(); }
	bool is_open() const{ return(out.is_open()); }

public:
	std::ofstream out;
};


}

#endif /* TWK_WRITER_H_ */
