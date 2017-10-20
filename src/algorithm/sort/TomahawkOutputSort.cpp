#include <cassert>

#include "TomahawkOutputSort.h"

namespace Tomahawk{
namespace Algorithm{
namespace Output{

bool TomahawkOutputSorter::sort(const std::string& input, const std::string& destinationPrefix, U64 memory_limit){
	if(!this->reader.Open(input)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << input << "..." << std::endl;
		return false;
	}

	//
	std::vector<std::string> paths = Helpers::filePathBaseExtension(destinationPrefix);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == Tomahawk::Constants::OUTPUT_LD_SUFFIX.size() &&
	   strncasecmp(&paths[3][0], &Tomahawk::Constants::OUTPUT_LD_SUFFIX[0], Tomahawk::Constants::OUTPUT_LD_SUFFIX.size()) == 0)
		this-> baseName = paths[2];
	else this->baseName = paths[1];

	// Writing
	this->reader.setWriterType(0);
	this->reader.addLiteral("\n##tomahawk_partialSortCommand=" + Helpers::program_string());
	this->reader.header.controller.sorted = 0;
	this->reader.header.controller.expanded = this->reverse_entries ? 1 : 0;
	this->reader.header.controller.partial_sort = 1;
	this->reader.OpenWriter(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX);

	basic_writer_type toi_writer;
	toi_writer.open(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX);
	toi_header_type headIndex(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, this->reader.header.samples, this->reader.header.n_contig);
	headIndex.controller.sorted = 0;
	headIndex.controller.expanded = this->reverse_entries ? 1 : 0;
	headIndex.controller.partial_sort = 1;
	toi_writer.getNativeStream() << headIndex;
	// writer
	basic_writer_type& stream = *reinterpret_cast<basic_writer_type*>(this->reader.writer->getStream());

	if(memory_limit < 10e6){
		memory_limit = 10e6;
		std::cerr << Helpers::timestamp("SORT") << "Setting memory limit to 10 MB..." << std::endl;
	}

	if(this->reverse_entries){
		std::cerr << Helpers::timestamp("SORT") << "Expanding: Reversing entries..." << std::endl;
		memory_limit /= 2;
	} else
		std::cerr << Helpers::timestamp("SORT") << "Maintaining: Not reversing..." << std::endl;

	// Perform indexed sorting if possible
	if(this->reader.hasIndex){
		return(this->__sortIndexed(toi_writer, input, memory_limit));
	}

	std::cerr << Helpers::timestamp("SORT") << "No index found..." << std::endl;

	bool trigger_break = false;
	U32 totempole_blocks_written = 0;
	while(true){
		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG","SORT") << "Reading..." << std::endl;

		if(!this->reader.nextBlockUntil(memory_limit))
			trigger_break = true;

		if(this->reader.output_buffer.size() == 0){
			trigger_break = true;
			break;
		}

		assert((this->reader.output_buffer.size() % sizeof(entry_type)) == 0);

		totempole_entry totempole;
		if(this->reverse_entries)
			totempole.entries = 2*(this->reader.output_buffer.size() / sizeof(entry_type));
		else
			totempole.entries = this->reader.output_buffer.size() / sizeof(entry_type);

		if(this->reverse_entries){
			const entry_type* entry = nullptr;
			while(this->reader.nextVariantLimited(entry)){
				entry_type temp(entry);
				temp.swapDirection();
				this->reader.output_buffer.Add((char*)&temp, sizeof(entry_type));
			}
		}

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG","SORT") << "Sorting: " << Helpers::ToPrettyString(this->reader.output_buffer.size()/sizeof(entry_sort_type)) << " entries" << std::endl;

		std::sort(reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[0]),
				  reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[this->reader.output_buffer.size()]));

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG","SORT") << "Indexing..." << std::endl;

		totempole.byte_offset = stream.getNativeStream().tellp();
		totempole.uncompressed_size = this->reader.output_buffer.size();

		this->reader.writer->write(this->reader.output_buffer);
		totempole.byte_offset_end = stream.getNativeStream().tellp();
		toi_writer.getNativeStream() << totempole;

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG","SORT") << "Writing..." << std::endl;

		++totempole_blocks_written;

		if(trigger_break) break;
	}

	// Make sure TOI is flushed before re-opening and seeking
	toi_writer.flush();

	// TOI
	// Update blocks written
	std::fstream re(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX, std::ios::in | std::ios::out | std::ios::binary);
	if(!re.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to reopen index..." << std::endl;
		return false;
	}

	re.seekg(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH + sizeof(float) + sizeof(U64) + sizeof(U32));
	if(!re.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to seek in index..." << std::endl;
		return false;
	}

	re.write((char*)&totempole_blocks_written, sizeof(U32));
	if(!re.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to update counts in index..." << std::endl;
		return false;
	}
	re.close();

	toi_writer.close();

	this->reader.writer->flush();
	this->reader.writer->close();

	return true;
}

/*
 * If we have an index available we can speed up sorting
 * by partitioning the data into N bins
 */
bool TomahawkOutputSorter::__sortIndexed(basic_writer_type& toi_writer, const std::string& input, U64 memory_limit){
	std::cerr << Helpers::timestamp("SORT") << "Index found..." << std::endl;

	std::vector< totempole_entry > blocks;

	totempole_entry totempole;
	totempole.byte_offset = this->reader.toi_reader[0].byte_offset;

	U64 n_entries = 0;
	for(U32 i = totempole.uncompressed_size; i < this->reader.toi_reader.size(); ++i){
		if(totempole.uncompressed_size > memory_limit){
			totempole.byte_offset_end = this->reader.toi_reader[i].byte_offset;
			blocks.push_back(totempole);
			totempole.byte_offset = this->reader.toi_reader[i].byte_offset;
			totempole.entries = 0;
			totempole.uncompressed_size = 0;
		}
		totempole.entries += this->reader.toi_reader[i].entries;
		totempole.uncompressed_size += this->reader.toi_reader[i].uncompressed_size;
		n_entries += totempole.entries;
	}

	// Have to add final
	if(blocks.size() == 0 || (totempole.byte_offset != blocks.back().byte_offset)){
		totempole.byte_offset_end = this->reader.toi_reader[this->reader.toi_reader.size() - 1].byte_offset_end;
		totempole.entries = 0;
		blocks.push_back(totempole);
	}

	if(totempole.entries != 0)
		blocks.push_back(totempole);

	// Split workload into different threads
	// Each thread get approximately 1/threads amount of work
	std::vector< totempole_entry > thread_workload(this->n_threads);
	const U64 limit_thread = (U64)((double)blocks.size()/this->n_threads);
	U32 current_thread_target = 0;
	U64 n_blocks_loaded = 0;
	totempole.reset();
	totempole.byte_offset = blocks[0].byte_offset;

	for(U32 i = 0; i < blocks.size(); ++i){
		if(n_blocks_loaded >= limit_thread && current_thread_target + 1 != this->n_threads){
			n_blocks_loaded = 0;
			totempole.byte_offset_end = blocks[i].byte_offset_end;
			thread_workload[current_thread_target] = totempole;
			totempole.byte_offset = blocks[i].byte_offset_end;
			++current_thread_target;
		}
		++n_blocks_loaded;
	}

	// Add final
	if(totempole.byte_offset != blocks.back().byte_offset_end){
		totempole.byte_offset_end = blocks.back().byte_offset_end;
		++current_thread_target;
		thread_workload[current_thread_target] = totempole;
	}

	std::cerr << Helpers::timestamp("SORT") << "Spawning " << current_thread_target << " threads..." << std::endl;
	std::thread** slaves = new std::thread*[current_thread_target];
	slave_sorter** instances = new slave_sorter*[current_thread_target];
	for(U32 i = 0; i < current_thread_target; ++i){
		instances[i] = new slave_sorter(this->reader.writer, toi_writer, memory_limit);
		if(!instances[i]->open(input)){
			std::cerr << Helpers::timestamp("ERROR", "SORT") << "Failed to reopen file..." << std::endl;
			exit(1);
		}

		// Trigger reverse if applicable
		instances[i]->reverseEntries(this->reverse_entries);
	}

	for(U32 i = 0; i < current_thread_target; ++i)
		slaves[i] = instances[i]->start(thread_workload[i]);

	for(U32 i = 0; i < current_thread_target; ++i)
		slaves[i]->join();

	U32 totempole_blocks_written = 0;
	for(U32 i = 0; i < current_thread_target; ++i)
		totempole_blocks_written += instances[i]->getBlocksWritten();

	// TOI
	// Update blocks written
	// Make sure TOI is flushed before re-opening and seeking
	toi_writer.flush();

	std::fstream re(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX, std::ios::in | std::ios::out | std::ios::binary);
	if(!re.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to reopen index..." << std::endl;
		return false;
	}

	re.seekg(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH + sizeof(float) + sizeof(U64) + sizeof(U32));
	if(!re.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to seek in index..." << std::endl;
		return false;
	}

	re.write((char*)&totempole_blocks_written, sizeof(U32));
	if(!re.good()){
		std::cerr << Helpers::timestamp("ERROR", "TWO") << "Failed to update counts in index..." << std::endl;
		return false;
	}
	re.close();

	toi_writer.close();

	this->reader.writer->flush();
	this->reader.writer->close();

	// Cleanup
	for(U32 i = 0; i < current_thread_target; ++i){
		delete instances[i];
		slaves[i] = nullptr;
	}
	delete [] instances;
	delete [] slaves;

	return true;
}

bool TomahawkOutputSorter::sortMerge(const std::string& inputFile, const std::string& destinationPrefix, const U32 block_size){
	if(!this->reader.Open(inputFile)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << inputFile << "..." << std::endl;
		return false;
	}

	if(!this->reader.hasIndex){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "File is not indexed!" << std::endl;
		return false;
	}

	this->reader.addLiteral("\n##tomahawk_mergeSortCommand=" + Helpers::program_string());

	toi_header_type& toi_header = this->reader.toi_reader.getHeader();
	if(toi_header.controller.sorted == true){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "File is already sorted!" << std::endl;
		return false;
	}

	if(toi_header.controller.partial_sort == false){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "File is not partially sorted!" << std::endl;
		return false;
	}

	toi_header.controller.partial_sort = false;
	toi_header.controller.sorted = true;
	this->reader.header.controller.partial_sort = false;
	this->reader.header.controller.sorted = true;
	writer_type writer(this->reader.contigs, &this->reader.header, toi_header);

	if(!writer.open(destinationPrefix)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open!" << std::endl;
		return false;
	}
	writer.writeHeader(this->reader.literals);

	const U32 n_toi_entries = this->reader.toi_reader.size();
	std::ifstream* streams = new std::ifstream[n_toi_entries];
	tgzf_iterator** iterators = new tgzf_iterator*[n_toi_entries];

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "SORT") << "Opening " << n_toi_entries << " file handles...";

	for(U32 i = 0; i < n_toi_entries; ++i){
		streams[i].open(inputFile);
		streams[i].seekg(this->reader.toi_reader[i].byte_offset);
		iterators[i] = new tgzf_iterator(streams[i], 65536, this->reader.toi_reader[i].byte_offset, this->reader.toi_reader[i].byte_offset_end);
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
		outQueue.push( queue_entry(e, i, IO::Support::TomahawkOutputEntryCompFuncConst) );
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
				outQueue.push( queue_entry(e, id, IO::Support::TomahawkOutputEntryCompFuncConst) );
				break;
			}
			writer << *e;
		}
	}

	writer.flush();
	if(!writer.finalize(toi_header.controller.expanded)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to finalize index..." << std::endl;
		return false;
	}

	writer.close();

	// Temp
	//index_type& index = writer.getIndex();
	//std::cerr << index << std::endl;

	// Cleanup
	for(U32 i = 0; i < n_toi_entries; ++i)
		delete iterators[i];

	delete [] iterators;
	delete [] streams;

	return true;
}

}
}
}
