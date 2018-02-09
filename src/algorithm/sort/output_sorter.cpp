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
	if(!this->reader.Open(input)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << input << "..." << std::endl;
		return false;
	}

	std::cerr << this->reader.filesize << "->" << this->reader.filesize / 8 << std::endl;
	const U64 n_variants_chunk = this->reader.filesize / this->n_threads;
	size_t n_total = 0;
	for(U32 i = 0; i < 8; ++i){
		OutputContainer o = this->reader.getContainerBytes(n_variants_chunk);
		std::cerr << "Size: " << o.size() << std::endl;
		if(o.size() == 0)
			continue;

		n_total += o.size();

		std::sort(&o.front(), &o.back(), &entry_type::sortAscending);

		const entry_type* prev = &o[0];
		std::cout << o[0] << '\n';
		for(size_t j = 1; j < o.size(); ++j){
			//std::cout.write((char*)&o[j], sizeof(entry_type));
			std::cout << o[j] << '\n';
			if(*prev >= o[j]){
				std::cerr << j-1 << ',' << j << std::endl;
				std::cerr << *prev << std::endl;
				std::cerr << o[j] << std::endl;
				exit(1);
			}
			prev = &o[j];
		}

		if(this->reverse_entries){
			// If we want to reverse the entries
			// Sketch:
			// 1: Swap A <> B
			// 2: Sort again
			for(size_t j = 1; j < o.size(); ++j)
				o[j].swapDirection();

			std::sort(&o.front(), &o.back(), &entry_type::sortAscending);

			prev = &o[0];
			std::cout << o[0] << '\n';
			for(size_t j = 1; j < o.size(); ++j){
				//std::cout.write((char*)&o[j], sizeof(entry_type));
				std::cout << o[j] << '\n';
				if(*prev >= o[j]){
					std::cerr << j-1 << ',' << j << std::endl;
					std::cerr << *prev << std::endl;
					std::cerr << o[j] << std::endl;
					exit(1);
				}
				prev = &o[j];
			}

		}
	}
	std::cerr << "total: " << n_total << std::endl;

	/*
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
	this->reader.OpenWriter(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX);

	// Open writer
	basic_writer_type toi_writer;
	toi_writer.open(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX);
	toi_header_type headIndex(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, this->reader.header.samples, this->reader.header.n_contig);
	headIndex.controller.sorted = 0;
	headIndex.controller.expanded = this->reverse_entries ? 1 : 0;
	headIndex.controller.partial_sort = 1;
	toi_writer.getNativeStream() << headIndex;
	// writer
	basic_writer_type& stream = *reinterpret_cast<basic_writer_type*>(this->reader.writer->getStream());

	if(memory_limit < 100e6){
		memory_limit = 100e6;
		std::cerr << Helpers::timestamp("SORT") << "Setting memory limit to 100 MB..." << std::endl;
	}

	if(this->reverse_entries){
		std::cerr << Helpers::timestamp("SORT") << "Reversing entries..." << std::endl;
	} else
		std::cerr << Helpers::timestamp("SORT") << "Not reversing..." << std::endl;

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

		if(this->reader.data_buffer.size() == 0){
			trigger_break = true;
			break;
		}

		assert((this->reader.data_buffer.size() % sizeof(entry_type)) == 0);

		totempole_entry totempole;
		if(this->reverse_entries)
			totempole.n_entries = 2*(this->reader.data_buffer.size() / sizeof(entry_type));
		else
			totempole.n_entries = this->reader.data_buffer.size() / sizeof(entry_type);

		if(this->reverse_entries){
			const entry_type* entry = nullptr;
			while(this->reader.nextVariantLimited(entry)){
				entry_type temp(entry);
				temp.swapDirection();
				this->reader.data_buffer.Add((char*)&temp, sizeof(entry_type));
			}
		}

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG","SORT") << "Sorting: " << Helpers::ToPrettyString(this->reader.data_buffer.size()/sizeof(entry_sort_type)) << " entries" << std::endl;

		std::sort(reinterpret_cast<entry_sort_type*>(this->reader.data_buffer.data()),
				  reinterpret_cast<entry_sort_type*>(&this->reader.data_buffer[this->reader.data_buffer.size()]));

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG","SORT") << "Indexing..." << std::endl;

		totempole.byte_offset = stream.getNativeStream().tellp();
		totempole.uncompressed_size = this->reader.data_buffer.size();

		this->reader.writer->write(this->reader.data_buffer);
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
	*/

	return true;
}

bool OutputSorter::__sortIndexed(basic_writer_type& toi_writer, const std::string& input, U64 memory_limit){
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
			totempole.n_entries = 0;
			totempole.uncompressed_size = 0;
		}
		totempole.n_entries += this->reader.toi_reader[i].n_entries;
		totempole.uncompressed_size += this->reader.toi_reader[i].uncompressed_size;
		n_entries += totempole.n_entries;
	}

	// Have to add final
	if(totempole.byte_offset != blocks.back().byte_offset){
		totempole.byte_offset_end = this->reader.toi_reader[this->reader.toi_reader.size() - 1].byte_offset_end;
		blocks.push_back(totempole);
	}

	if(totempole.n_entries != 0)
		blocks.push_back(totempole);

	// Todo: if n_threads > blocks.size()
	// set n_threads to block.size() and give each thread 1 block

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
			totempole.byte_offset_end = blocks[i].byte_offset;
			thread_workload[current_thread_target] = totempole;
			totempole.byte_offset = blocks[i].byte_offset;
			++current_thread_target;
		}
		++n_blocks_loaded;
	}

	// Add final
	if(current_thread_target != this->n_threads){
		totempole.byte_offset_end = blocks.back().byte_offset_end;
		thread_workload[this->n_threads-1] = totempole;
	}

	std::cerr << Helpers::timestamp("SORT") << "Spawning " << this->n_threads << " threads..." << std::endl;
	std::thread** slaves = new std::thread*[this->n_threads];
	slave_sorter** instances = new slave_sorter*[this->n_threads];
	for(U32 i = 0; i < this->n_threads; ++i){
		instances[i] = new slave_sorter(this->reader.writer, toi_writer, memory_limit);
		if(!instances[i]->open(input)){
			std::cerr << Helpers::timestamp("ERROR", "SORT") << "Failed to reopen file..." << std::endl;
			exit(1);
		}

		// Trigger reverse if applicable
		instances[i]->reverseEntries(this->reverse_entries);
	}

	for(U32 i = 0; i < this->n_threads; ++i)
		slaves[i] = instances[i]->start(thread_workload[i]);

	for(U32 i = 0; i < this->n_threads; ++i)
		slaves[i]->join();

	U32 totempole_blocks_written = 0;
	for(U32 i = 0; i < this->n_threads; ++i)
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
	for(U32 i = 0; i < this->n_threads; ++i){
		delete instances[i];
		slaves[i] = nullptr;
	}
	delete [] instances;
	delete [] slaves;

	return true;
}

bool OutputSorter::sortMerge(const std::string& inputFile, const std::string& destinationPrefix, const U32 block_size){
	/*
	if(!this->reader.Open(inputFile)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << inputFile << "..." << std::endl;
		return false;
	}

	if(!this->reader.hasIndex){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "File is not indexed!" << std::endl;
		return false;
	}

	this->reader.addLiteral("\n##tomahawk_mergeSortCommand=" + Helpers::program_string());

	toi_header_type toi_header = this->reader.toi_reader.getHeader();
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
	8
	*/

	return true;
}

}
}
