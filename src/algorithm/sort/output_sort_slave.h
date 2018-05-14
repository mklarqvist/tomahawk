#ifndef TOMAHAWKOUTPUTSORTSLAVE_H_
#define TOMAHAWKOUTPUTSORTSLAVE_H_

#include <cassert>

#include "io/output_writer.h"

namespace tomahawk{
namespace algorithm{

/**<
 * Worker slave for partially sorting a (large) `two` file
 */
class OutputSortSlave {
private:
	typedef OutputSortSlave       self_type;
	typedef io::OutputEntry       entry_type;
	typedef io::OutputWriter      writer_type;
	typedef TomahawkOutputReader  reader_type;
	typedef io::TGZFController    tgzf_controller_type;
	typedef io::BasicBuffer       buffer_type;
	typedef io::TGZFEntryIterator<entry_type> tgzf_iterator;

public:
	OutputSortSlave(reader_type& reader, writer_type& writer, const std::pair<U64, U64>& workload, const U32 memory_limit) :
		workload_(workload),
		n_memory_limit_(memory_limit),
		reader_(reader),
		writer_(writer)
	{}

	~OutputSortSlave(){
		this->inflate_buffer_.deleteAll();
		this->data_.deleteAll();
	}

	bool open(const std::string& input){
		this->stream_.open(input, std::ios::binary | std::ios::in);
		if(this->stream_.good() == false){
			std::cerr << helpers::timestamp("ERROR","SORT") << "Failed to open: " << input << "..." << std::endl;
			return false;
		}

		return true;
	}

	std::thread* start(void){
		this->thread_ = std::thread(&self_type::sort, this);
		return(&this->thread_);
	}

	inline const writer_type& getWriter(void) const{ return(this->writer_); }

private:
	bool sort(void){
		if(!this->stream_.seekg(this->workload_.first)){
			std::cerr << helpers::timestamp("ERROR","SORT") << "Failed to seek to block!" << std::endl;
			exit(1); // exit instead of return because of detached threads
		}

		// iterator
		const U64 n_entries_limit = this->n_memory_limit_ / sizeof(entry_type);
		tgzf_iterator it(this->stream_, 65536, this->workload_.first, this->workload_.second);
		bool finished_ = false;

		if(!this->stream_.good()){
			std::cerr << helpers::timestamp("ERROR","SORT") << "Stream is bad!" << std::endl;
			exit(1); // exit instead of return because of detached threads
		}

		while(true){
			OutputContainer container(this->n_memory_limit_ / sizeof(entry_type) + 1024);

			const entry_type* e = nullptr;
			for(U32 i = 0; i < n_entries_limit; ++i){
				if(!it.nextEntry(e)){
					finished_ = true;
					break;
				}
				container += *e;
			}
			std::sort(&container.front(), &container.back() + 1);

			const entry_type* prev = &container[0];
			for(size_t j = 1; j < container.size(); ++j){
				if(*prev >= container[j]){
					std::cerr << helpers::timestamp("ERROR","SORT");
					std::cerr << j-1 << ',' << j << std::endl;
					std::cerr << *prev << std::endl;
					std::cerr << container[j] << std::endl;
					exit(1);
				}
				prev = &container[j];
			}

			this->writer_ << container;

			if(finished_)
				break;
		}

		return true;
	}

private:
	std::ifstream        stream_;
	std::pair<U64, U64>  workload_;
	const U32            n_memory_limit_;
	const reader_type&   reader_;
	writer_type          writer_;
	std::thread          thread_;
	buffer_type          inflate_buffer_;
	buffer_type          data_;
	tgzf_controller_type compression_manager_;
};


}
}

#endif /* TOMAHAWKOUTPUTSORTSLAVE_H_ */
