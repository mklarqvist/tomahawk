#ifndef TWK_WRITER_H_
#define TWK_WRITER_H_

#include <ostream>
#include "index.h"
#include "core.h"
#include "two_reader.h"
#include "spinlock.h"

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
		//std::cerr << "wrote=" << b_unc << "->" << b_comp << " " << (float)b_unc/b_comp << " -> " << obuf.size() << std::endl;

		spinlock.unlock();
		obuf.reset();
	}

	void Add(const uint32_t b_unc, const uint32_t b_comp, twk_buffer_t& obuf, IndexEntryOutput& entry){
		spinlock.lock();
		stream.flush();

		entry.foff = stream.tellp();

		uint8_t marker = 1; // non-termination marker
		SerializePrimitive(marker, stream);
		SerializePrimitive(b_unc,  stream); // uncompressed size
		SerializePrimitive(b_comp, stream); // compressed size
		stream.write(obuf.data(), obuf.size()); // actual data
		//std::cerr << "wrote=" << b_unc << "->" << b_comp << " " << (float)b_unc/b_comp << " -> " << obuf.size() << std::endl;
		stream.flush();
		entry.fend = stream.tellp();

		spinlock.unlock();
		obuf.reset();
	}

	static std::string RandomSuffix() {
		 std::string name1 = std::tmpnam(nullptr);
		 std::vector<std::string> ret = utility::split(name1, '/');
		 return(ret.back());
	}

	static std::string GetExtension(const std::string ret) {
		return(utility::ExtensionName(ret));
	}

	static std::string GetBasePath(const std::string ret) {
		return(utility::BasePath(ret));
	}

	static std::string GetBaseName(const std::string ret) {
		return(utility::BaseName(ret));
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


struct twk_two_writer_t : public twk_writer_t {
	twk_two_writer_t() : mode('u'), c_level(1), n_blk_lim(10000), oblock(10000){
		buf = std::cout.rdbuf();
		stream.basic_ios<char>::rdbuf(buf);
	}

	inline void SetCompressionLevel(const int32_t level){ c_level = level; }

	bool Open(const std::string& file){
		if(file.size() == 0 || (file.size() == 1 && file[0] == '-')) return _OpenStream();
		else return _OpenFile(file);
	}

	bool _OpenFile(const std::string& file){
		out.open(file,std::ios::out | std::ios::binary);
		if(!out.good()){
			std::cerr << "failed to open" << std::endl;
			return false;
		}
		buf = out.rdbuf();
		stream.rdbuf(buf);
		return true;
	}

	bool _OpenStream(){
		buf = std::cout.rdbuf();
		stream.basic_ios<char>::rdbuf(buf);
		return true;
	}

	inline void close(){ out.close(); }
	inline bool is_open() const{ return true; }

	bool WriteHeader(two_reader& reader){
		if(mode == 'b') return(WriteHeaderBinary(reader));
		else if(mode == 'u') return(WriteHeaderReadable(reader));
		return false;
	}

	bool WriteHeaderBinary(two_reader& reader) {
		if(stream.good() == false)
			return false;

		obuf.reset();
		stream.write(tomahawk::TOMAHAWK_LD_MAGIC_HEADER.data(), tomahawk::TOMAHAWK_LD_MAGIC_HEADER_LENGTH);
		tomahawk::twk_buffer_t buf(256000);
		buf << reader.hdr;
		if(zcodec.Compress(buf, obuf, c_level) == false){
			std::cerr << "failed to compress" << std::endl;
			return false;
		}

		stream.write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		stream.write(obuf.data(),obuf.size());
		stream.flush();
		buf.reset();

		return(stream.good());
	}

	bool WriteHeaderReadable(two_reader& reader) {
		if(stream.good() == false)
			return false;

		obuf.reset();
		stream << reader.hdr.literals_;
		stream << "flags\tridA\tposA\tridB\tposB\tHOMHOM\tHOMALT\tALTHOM\tALTALT\tD\tDprime\tR\tR2\tP\tChiSqFisher\tChiSqModel" << std::endl;
		stream.flush();

		return(stream.good());
	}

	bool WriteFinal(){
		//std::cerr << "in write final=" << stream.good() << " and mode=" << mode << std::endl;
		if(stream.good() == false) return false;
		if(WriteBlock() == false) return false;

		obuf.reset();
		tomahawk::twk_buffer_t buf(256000);

		buf << oindex;
		//std::cerr << "index buf size =" << buf.size() << std::endl;
		if(zcodec.Compress(buf, obuf, c_level) == false){
			std::cerr << "failed compression" << std::endl;
			return false;
		}
		//std::cerr << buf.size() << "->" << obuf.size() << " -> " << (float)buf.size()/obuf.size() << std::endl;

		// temp write out index
		//for(int i = 0; i < oindex.n; ++i){
		//	std::cerr << i << "/" << oindex.n << " " << oindex.ent[i].rid << ":" << oindex.ent[i].minpos << "-" << oindex.ent[i].maxpos << " offset=" << oindex.ent[i].foff << "-" << oindex.ent[i].fend << " brid=" << oindex.ent[i].ridB << std::endl;
		//}

		//for(int i = 0; i < oindex.m_ent; ++i){
		//	std::cerr << "meta=" << i << "/" << oindex.m_ent << " " << oindex.ent_meta[i].rid << ":" << oindex.ent_meta[i].minpos << "-" << oindex.ent_meta[i].maxpos << " offset=" << oindex.ent_meta[i].foff << "-" << oindex.ent_meta[i].fend << std::endl;
		//}

		const uint64_t offset_start_index = stream.tellp();
		uint8_t marker = 0;
		stream.write(reinterpret_cast<const char*>(&marker),sizeof(uint8_t));
		stream.write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		stream.write(obuf.data(),obuf.size());
		stream.write(reinterpret_cast<const char*>(&offset_start_index),sizeof(uint64_t));
		stream.write(tomahawk::TOMAHAWK_FILE_EOF.data(), tomahawk::TOMAHAWK_FILE_EOF_LENGTH);
		stream.flush();
		return(stream.good());
	}

	bool Add(const twk1_two_t& rec){
		if(oblock.n == n_blk_lim){
			if(WriteBlock() == false)
				return false;
		}

		oblock += rec;
		return true;
	}

	bool WriteBlock(){
		if(mode == 'b') return(WriteBlockCompressedTWO());
		else if(mode == 'u') return(WriteBlockUncompressedLD());

		return false;
	}

	bool WriteBlockUncompressedLD(){
		if(oblock.n){
			for(int i = 0; i < oblock.n; ++i){
				oblock.rcds[i].PrintLD(std::cout);
			}
			oblock.reset();
			std::cout.flush();
		}
		return true;
	}

	bool WriteBlockCompressedTWO(){
		// Todo: if uncompressed data then write as is
		if(oblock.n){
			//twk_oblock_two_t b;
			ubuf << oblock;

			if(zcodec.Compress(ubuf, obuf, c_level) == false){
				std::cerr << "failed compression" << std::endl;
				return false;
			}

			// parent add
			ioentry.foff = stream.tellp();
			twk_writer_t::Add(ubuf.size(), obuf.size(), obuf);

			ioentry.n     = oblock.n;
			ioentry.b_unc = ubuf.size();
			ioentry.b_cmp = obuf.size();

			if(oindex.state == TWK_IDX_SORTED){ // if index is sorted
				uint32_t ridb = oblock.rcds[0].ridB;
				for(int i = 1; i < oblock.n; ++i){ // check if ridb is uniform
					if(oblock.rcds[i].ridB != ridb){
						ridb = -1;
						break;
					}
				}
				ioentry.rid    = oblock.rcds[0].ridA;
				ioentry.ridB   = ridb;
				ioentry.minpos = oblock.rcds[0].Apos;
				ioentry.maxpos = oblock.rcds[oblock.n-1].Apos;
			} else {
				ioentry.rid    = -1;
				ioentry.ridB   = -1;
				ioentry.minpos = 0;
				ioentry.maxpos = 0;
			}
			ioentry.fend = stream.tellp();

			oindex += ioentry;
			if(oindex.state == TWK_IDX_SORTED){
				oindex.ent_meta[ioentry.rid] += ioentry;
			}

			ioentry.clear();
			obuf.reset(); ubuf.reset();
			oblock.reset();
		}
		return true;
	}

	char mode;
	int32_t c_level;
	uint32_t n_blk_lim; // flush block limit
	twk1_two_block_t oblock;
	IndexEntryOutput ioentry;
	IndexOutput oindex;
	ZSTDCodec zcodec;
	twk_buffer_t ubuf, obuf;
	std::ofstream out;
};

}

#endif /* TWK_WRITER_H_ */
