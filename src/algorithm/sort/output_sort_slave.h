#ifndef TOMAHAWKOUTPUTSORTSLAVE_H_
#define TOMAHAWKOUTPUTSORTSLAVE_H_

#include <cassert>

#include "../../io/output_writer.h"

namespace Tomahawk{
namespace Algorithm{

/**<
 * Worker slave for partially sorting a (large) `two` file
 */
class OutputSortSlave {
private:
	typedef OutputSortSlave       self_type;
	typedef IO::OutputEntry       entry_type;
	typedef IO::OutputWriter      writer_type;
	typedef TomahawkOutputReader  reader_type;
	typedef IO::TGZFController    tgzf_controller_type;
	typedef IO::BasicBuffer       buffer_type;

public:
	OutputSortSlave(reader_type& reader, writer_type& writer, const std::pair<U32, U32>& workload, const U32 memory_limit) :
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
			std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << input << "..." << std::endl;
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
		if(!this->reader_.seekBlock(this->stream_, this->workload_.first)){
			std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to seek to block!" << std::endl;
			exit(1); // exit instead of return because of detached threads
		}

		const U32 n_blocks = this->workload_.second - this->workload_.first;
		for(U32 i = 0; i < n_blocks; ++i){
			if(!this->reader_.parseBlock(this->stream_, this->inflate_buffer_, this->data_, this->compression_manager_, false)){
				std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to get block!" << std::endl;
				exit(1); // exit instead of return because of detached threads
			}
		}

		// Load a chunk of data
		OutputContainer o(this->data_.data(), this->data_.size());
		if(o.size() == 0)
			return true;

		std::sort(&o.front(), &o.back());

		const entry_type* prev = &o[0];
		for(size_t j = 1; j < o.size(); ++j){
			if(*prev >= o[j]){
				std::cerr << j-1 << ',' << j << std::endl;
				std::cerr << *prev << std::endl;
				std::cerr << o[j] << std::endl;
				exit(1);
			}
			prev = &o[j];
		}

		this->writer_ << o;

		return true;
	}

private:
	std::ifstream        stream_;
	std::pair<U32, U32>  workload_;
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
