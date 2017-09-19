#ifndef TOMAHAWKOUTPUTSORTSLAVE_H_
#define TOMAHAWKOUTPUTSORTSLAVE_H_

#include <cassert>

namespace Tomahawk{
namespace Algorithm{
namespace Output{

class TomahawkOutputSortSlave{
	typedef TomahawkOutputSortSlave self_type;
	typedef IO::TomahawkOutputEntry entry_type;
	typedef IO::TomahawkOutputEntrySort entry_sort_type;
	typedef Totempole::TotempoleOutputEntry totempole_entry;
	typedef IO::WriterFile writer_type;
	typedef IO::TomahawkOutputReader two_reader_type;
	typedef Tomahawk::IO::TomahawkOutputWriterInterface two_writer_interface;
	typedef Tomahawk::IO::TomahawkOutputWriter two_writer_type;
	typedef IO::TGZFController tgzf_controller_type;

public:
	TomahawkOutputSortSlave(two_writer_interface* writer, writer_type& toi_writer, const U32 memory_limit) :
		memory_limit(memory_limit),
		writer(reinterpret_cast<two_writer_type*>(writer)),
		toi_writer(toi_writer)
	{}
	~TomahawkOutputSortSlave(){}

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


			if(!this->reader.nextBlockUntil(memory_limit))
				trigger_break = true;

			if(this->reader.output_buffer.size() == 0){
				trigger_break = true;
				break;
			}

			assert((this->reader.output_buffer.size() % sizeof(entry_type)) == 0);

			std::sort(reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[0]),
					  reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[this->reader.output_buffer.size()]));

			const entry_type* entry;
			if(!this->reader.nextVariantLimited(entry)){
				std::cerr << Helpers::timestamp("ERROR","SORT") << "No data!" << std::endl;
				return false;
			}

			totempole.reset();
			totempole.entries = 1;
			totempole.uncompressed_size = this->reader.output_buffer.size();

			const entry_type* prev;
			std::swap(prev, entry);

			while(this->reader.nextVariantLimited(entry)){
				++totempole.entries;
				std::swap(prev, entry);
			}

			this->controller.Clear();
			this->controller.Deflate(this->reader.output_buffer);

			this->writer->getLock()->lock();
			totempole.byte_offset = stream.getNativeStream().tellp();
			this->writer->getStream()->writeNoLock(this->controller.buffer.data, this->controller.buffer.pointer);
			totempole.byte_offset_end = stream.getNativeStream().tellp();
			toi_writer.getNativeStream() << totempole;
			this->writer->getLock()->unlock();

			if(trigger_break) break;
		}

		return true;
	}

private:
	const U32 memory_limit;
	two_writer_type* writer;
	writer_type& toi_writer;
	two_reader_type reader;
	std::thread thread;
	tgzf_controller_type controller;
};


}
}
}

#endif /* TOMAHAWKOUTPUTSORTSLAVE_H_ */
