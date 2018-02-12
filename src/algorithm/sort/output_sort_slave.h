#ifndef TOMAHAWKOUTPUTSORTSLAVE_H_
#define TOMAHAWKOUTPUTSORTSLAVE_H_

#include <cassert>

namespace Tomahawk{
namespace Algorithm{

/*
class OutputSortSlave{
private:
	typedef OutputSortSlave                   self_type;
	typedef IO::OutputEntry                   entry_type;
	typedef IO::WriterFile                    writer_type;
	typedef IO::TomahawkOutputReader          two_reader_type;
	typedef IO::TGZFController                tgzf_controller_type;

public:

	OutputSortSlave(two_writer_interface* writer, writer_type& toi_writer, const U32 memory_limit) :
		memory_limit(memory_limit),
		blocks_written(0),
		writer(reinterpret_cast<two_writer_type*>(writer)),
		toi_writer(toi_writer),
		reverse_entries(true)
	{}
	~OutputSortSlave(){}

	inline void reverseEntries(const bool yes = true){ this->reverse_entries = yes; }

	bool open(const std::string& input){
		if(!this->reader.Open(input)){
			std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << input << "..." << std::endl;
			return false;
		}

		return true;
	}

	std::thread* start(const totempole_entry& workload){
		this->thread = std::thread(&self_type::sort, this, workload);
		return(&this->thread);
	}

	inline const U32& getBlocksWritten(void) const{ return(this->blocks_written); }


private:

	bool sort(const totempole_entry& workload){

		bool trigger_break = false;
		writer_type& stream = *reinterpret_cast<writer_type*>(this->writer->getStream());
		this->reader.stream.seekg(workload.byte_offset);
		if(!this->reader.stream.good()){
			std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to seek in file..." << std::endl;
			return false;
		}

		totempole_entry totempole;
		while(true){
			if(this->reader.stream.tellg() == workload.byte_offset_end)
				break;

			if(!this->reader.nextBlockUntil(this->memory_limit, workload.byte_offset_end))
				trigger_break = true;


			if(this->reader.data_buffer.size() == 0){
				trigger_break = true;
				break;
			}

			assert((this->reader.data_buffer.size() % sizeof(entry_type)) == 0);

			if(this->reverse_entries){
				const entry_type* entry = nullptr;
				totempole.n_entries += 2*((this->reader.data_buffer.size() % sizeof(entry_type)));
				while(this->reader.nextVariantLimited(entry)){
					// Flip cA,pA with cB,pB
					entry_type temp(entry);
					temp.swapDirection();
					this->reader.data_buffer.Add((char*)&temp, sizeof(entry_type));
				}
			} else {
				// Do not reverse
				totempole.n_entries = (this->reader.data_buffer.size() % sizeof(entry_type));
			}

			std::sort(reinterpret_cast<entry_sort_type*>(this->reader.data_buffer.data()),
					  reinterpret_cast<entry_sort_type*>(&this->reader.data_buffer[this->reader.data_buffer.size()]));

			totempole.reset();
			totempole.uncompressed_size = this->reader.data_buffer.size();

			this->controller.Clear();
			this->controller.Deflate(this->reader.data_buffer);

			this->writer->getLock()->lock();
			totempole.byte_offset = stream.getNativeStream().tellp();
			this->writer->getStream()->writeNoLock(this->controller.buffer.data(), this->controller.buffer.size());
			++this->blocks_written;
			totempole.byte_offset_end = stream.getNativeStream().tellp();
			toi_writer.getNativeStream() << totempole;
			this->writer->getLock()->unlock();

			if(trigger_break) break;
		}

		return true;
	}


private:
	const U32 memory_limit;
	U32 blocks_written;
	two_writer_type* writer;
	writer_type& toi_writer;
	two_reader_type reader;
	std::thread thread;
	tgzf_controller_type controller;
	bool reverse_entries;
};
*/

}
}

#endif /* TOMAHAWKOUTPUTSORTSLAVE_H_ */
