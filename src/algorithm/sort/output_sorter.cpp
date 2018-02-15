#include <cassert>

#include "output_sorter.h"

namespace Tomahawk{
namespace Algorithm{

// Algorithmic sketch:
// 1: Load data into balanced chunks of memory_limit bytes
// 2: Load partitioned data into containers
// 3: Perform sort
// 4: Perform merge if desired
bool OutputSorter::sort(const std::string& input, const std::string& destinationPrefix, U64 memory_limit){
	if(!this->reader.open(input)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << input << "..." << std::endl;
		return false;
	}

	if(this->reader.getIndex().getController().isSorted){
		std::cerr << Helpers::timestamp("LOG","SORT") << "File is already sorted..." << std::endl;
		return true;
	}

	//std::cerr << this->reader.getIndex().totalBytes() << "->" << this->reader.getIndex().totalBytes() / this->n_threads << std::endl;
	const U64 n_variants_chunk = this->reader.getIndex().totalBytes() / this->n_threads;

	std::pair<U32, U32>* thread_distribution = new std::pair<U32, U32>[this->n_threads];
	size_t i = 0;
	for(U32 t = 0; t < this->n_threads; ++t){
		U64 partition_size = 0;
		const size_t from = i;
		for(; i < this->reader.getIndex().size(); ++i){
			partition_size += this->reader.getIndex().getContainer().at(i).sizeBytes();
			if(partition_size >= n_variants_chunk && t + 1 != this->n_threads){
				++i;
				break;
			}
		}
		thread_distribution[t].first = from;
		thread_distribution[t].second = i;

		//std::cerr << "t: " << from << "->" << i << " for " << partition_size << "/" << n_variants_chunk << std::endl;
		if(i == this->reader.getIndex().size()){
			//std::cerr << "ran out of data" << std::endl;
			break;
		}
	}

	// Append executed command to literals
	this->reader.getHeader().getLiterals() += "\n##tomahawk_sortCommand=" + Helpers::program_string();

	IO::OutputWriter writer;
	if(!writer.open(destinationPrefix)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << destinationPrefix << "..." << std::endl;
		return false;
	}
	writer.WriteHeaders(this->reader.getHeader());

	OutputSortSlave** slaves = new OutputSortSlave*[this->n_threads];
	std::thread** threads = new std::thread*[this->n_threads];

	for(U32 i = 0; i < this->n_threads; ++i){
		slaves[i] = new OutputSortSlave(this->reader, writer, thread_distribution[i], 1);
		if(!slaves[i]->open(input)){
			std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << input << "..." << std::endl;
			return false;
		}
		threads[i] = slaves[i]->start();
	}

	for(U32 i = 0; i < this->n_threads; ++i)
		threads[i]->join();

	for(U32 i = 0; i < this->n_threads; ++i)
		writer += slaves[i]->getWriter();

	writer.setSorted(false);
	writer.setPartialSorted(true);
	writer.flush();
	writer.WriteFinal();

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG") << "Output: " << Helpers::ToPrettyString(writer.sizeEntries()) << " entries into " << Helpers::ToPrettyString(writer.sizeBlocks()) << " blocks..." << std::endl;

	for(U32 i = 0; i < this->n_threads; ++i)
		delete slaves[i];

	delete [] slaves;
	delete [] threads;
	delete [] thread_distribution;

	return true;
}

bool OutputSorter::sortMerge(const std::string& inputFile, const std::string& destinationPrefix, const U32 block_size){
	if(!this->reader.open(inputFile)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << inputFile << "..." << std::endl;
		return false;
	}

	if(this->reader.getIndex().getController().isSorted){
		std::cerr << Helpers::timestamp("LOG","SORT") << "File is already sorted..." << std::endl;
		return true;
	}

	if(this->reader.getIndex().getController().isPartialSorted == false){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "File is not partially sorted..." << std::endl;
		return false;
	}

	// Append executed command to literals
	this->reader.getHeader().getLiterals() += "\n##tomahawk_mergeSortCommand=" + Helpers::program_string();

	IO::OutputWriter writer;
	if(!writer.open(destinationPrefix)){
		std::cerr << Helpers::timestamp("ERROR", "SORT") << "Failed to open: " << destinationPrefix << "..." << std::endl;
		return false;
	}
	writer.WriteHeaders(this->reader.getHeader());

	const U32 n_toi_entries = this->reader.getIndex().size();
	std::ifstream* streams = new std::ifstream[n_toi_entries];
	tgzf_iterator** iterators = new tgzf_iterator*[n_toi_entries];

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "SORT") << "Opening " << n_toi_entries << " file handles...";

	for(U32 i = 0; i < n_toi_entries; ++i){
		streams[i].open(inputFile);
		streams[i].seekg(this->reader.getIndex().getContainer()[i].byte_offset);
		iterators[i] = new tgzf_iterator(streams[i], 65536, this->reader.getIndex().getContainer()[i].byte_offset, this->reader.getIndex().getContainer()[i].byte_offset_end);
	}

	if(!SILENT)
		std::cerr << " Done!" << std::endl;

	// queue
	queue_type outQueue;

	//
	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "SORT") << "Merging..." << std::endl;

	// draw one from each
	const entry_type* e = nullptr;
	for(U32 i = 0; i < n_toi_entries; ++i){
		if(!iterators[i]->nextEntry(e)){
			std::cerr << Helpers::timestamp("ERROR", "SORT") << "Failed to get an entry..." << std::endl;
			return false;
		}
		outQueue.push( queue_entry(e, i, entry_type::sortAscending) );
	}

	if(outQueue.empty()){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "No data in queue..." << std::endl;
		return false;
	}

	// while queue is not empty
	while(outQueue.empty() == false){
		// peek at top entry in queue
		const U32 id = outQueue.top().streamID;
		writer << outQueue.top().data;

		// remove this record from the queue
		outQueue.pop();


		while(iterators[id]->nextEntry(e)){
			if(!(*e < outQueue.top().data)){
				outQueue.push( queue_entry(e, id, entry_type::sortAscending) );
				break;
			}
			writer << *e;
		}
	}

	writer.setPartialSorted(false);
	writer.setSorted(true);
	writer.flush();
	writer.WriteFinal();

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG") << "Output: " << Helpers::ToPrettyString(writer.sizeEntries()) << " entries into " << Helpers::ToPrettyString(writer.sizeBlocks()) << " blocks..." << std::endl;

	// Cleanup
	for(U32 i = 0; i < n_toi_entries; ++i)
		delete iterators[i];

	delete [] iterators;
	delete [] streams;

	return true;
}

}
}
