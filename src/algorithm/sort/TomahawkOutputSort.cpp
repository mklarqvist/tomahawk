#include <cassert>

#include "TomahawkOutputSort.h"

namespace Tomahawk{
namespace Algorithm{
namespace Output{

bool TomahawkOutputSorter::__sortIndexed(basic_writer_type& toi_writer, const std::string& input, const U32 memory_limit){
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

	if(totempole.entries != 0)
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
			std::cerr << Helpers::timestamp("ERROR", "SORT") << "Failed reopening file..." << std::endl;
			exit(1);
		}
	}

	for(U32 i = 0; i < this->n_threads; ++i){
		slaves[i] = instances[i]->start(thread_workload[i]);
	}

	for(U32 i = 0; i < this->n_threads; ++i)
		slaves[i]->join();

	toi_writer.flush();
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

bool TomahawkOutputSorter::sort(const std::string& input, const std::string& destinationPrefix, const U64 memory_limit){
	if(!this->reader.Open(input)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << input << "..." << std::endl;
		return false;
	}

	//
	std::vector<std::string> paths = Helpers::filePathBaseExtension(destinationPrefix);
	std::string basePath = paths[0];
	if(basePath.size() > 0)
		basePath += '/';

	std::string baseName;
	if(paths[3].size() == Tomahawk::Constants::OUTPUT_LD_SUFFIX.size() &&
	   strncasecmp(&paths[3][0], &Tomahawk::Constants::OUTPUT_LD_SUFFIX[0], Tomahawk::Constants::OUTPUT_LD_SUFFIX.size()) == 0)
		 baseName = paths[2];
	else baseName = paths[1];

	// Writing
	this->reader.setWriterType(0);
	this->reader.addLiteral("\n##tomahawk_partialSortCommand=" + Helpers::program_string());
	this->reader.OpenWriter(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX);
	basic_writer_type toi_writer;
	toi_writer.open(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX);
	IO::TomahawkOutputSortHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> headIndex(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, this->reader.header.samples, this->reader.header.n_contig);
	toi_writer.getNativeStream() << headIndex;
	// writer
	basic_writer_type& stream = *reinterpret_cast<basic_writer_type*>(this->reader.writer->getStream());

	//
	if(this->reader.hasIndex){
		return(this->__sortIndexed(toi_writer, input, memory_limit));
	}

	std::cerr << Helpers::timestamp("SORT") << "No index found..." << std::endl;

	bool trigger_break = false;
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

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG","SORT") << "Sorting: " << Helpers::ToPrettyString(this->reader.output_buffer.size()/sizeof(entry_sort_type)) << " entries" << std::endl;

		std::sort(reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[0]),
				  reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[this->reader.output_buffer.size()]));

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG","SORT") << "Indexing..." << std::endl;

		const entry_type* entry;
		if(!this->reader.nextVariantLimited(entry)){
			std::cerr << Helpers::timestamp("ERROR","SORT") << "No data!" << std::endl;
			return false;
		}
		totempole_entry totempole;
		totempole.entries = 1;
		totempole.byte_offset = stream.getNativeStream().tellp();
		totempole.uncompressed_size = this->reader.output_buffer.size();

		const entry_type* prev;
		std::swap(prev, entry);

		while(this->reader.nextVariantLimited(entry)){
			++totempole.entries;
			std::swap(prev, entry);
		}

		this->reader.writer->write(this->reader.output_buffer);
		totempole.byte_offset_end = stream.getNativeStream().tellp();
		toi_writer.getNativeStream() << totempole;

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG","SORT") << "Writing..." << std::endl;

		if(trigger_break) break;
	}

	toi_writer.flush();
	toi_writer.close();

	this->reader.writer->flush();
	this->reader.writer->close();

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

	writer_type writer(this->reader.contigs, &this->reader.header);
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
	const entry_type* e;
	for(U32 i = 0; i < n_toi_entries; ++i){
		if(!iterators[i]->nextEntry(e)){
			std::cerr << Helpers::timestamp("ERROR", "SORT") << "Failed to get an entry..." << std::endl;
			return false;
		}
		//entry_type hard_copy(e); // invoke copy ctor to avoid pointer errors when modifying internal buffer
		outQueue.push( queue_entry(e, i, IO::Support::TomahawkOutputEntryCompFuncConst) );
	}

	if(outQueue.empty()){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "No data in queue..." << std::endl;
		return false;
	}

	// Secondary TOI index
	/*
	U32 currentAID = outQueue.top().data.AcontigID;
	U32 currentAPos = outQueue.top().data.Aposition;
	U64 AIDSteps = 0;
	U64 APosSteps = 0;
	double AposStepsR = 0;
	U64 outputEntries = 0;
	*/
	//

	// while queue is not empty
	while(outQueue.empty() == false){
		// peek at top entry in queue
		const U32 id = outQueue.top().streamID;
		writer << outQueue.top().data;

		//std::cerr << outQueue.top().data.Aposition / (this->reader.contigs[outQueue.top().data.AcontigID].bases >> 9) << std::endl;
		//assert((outQueue.top().data.Aposition / (this->reader.contigs[outQueue.top().data.AcontigID].bases >> 9)) < 1024);

		//
		/*
		const entry_type& ent = outQueue.top().data;
		++AIDSteps;
		++APosSteps;
		AposStepsR += ent.R2;

		if(ent.Aposition != currentAPos || ent.AcontigID != currentAID){
			std::cout << "2switch: " << currentAID << '\t' << currentAPos << '\t' << APosSteps << '\t' << AposStepsR/APosSteps << '\n';
			currentAPos = ent.Aposition;
			APosSteps = 0;
			AposStepsR = 0;
			++outputEntries;
		}

		if(ent.AcontigID != currentAID){
			std::cerr << "1switch: " << currentAID << "->" << ent.AcontigID << '\t' << AIDSteps << std::endl;
			currentAID = ent.AcontigID;
			AIDSteps = 0;
			++outputEntries;
		}
		*/
		//

		// remove this record from the queue
		outQueue.pop();

		// Replace value from target stream
		//test
		//int c = 0;


		while(iterators[id]->nextEntry(e)){
			if(!(*e < outQueue.top().data)){
				outQueue.push( queue_entry(e, id, IO::Support::TomahawkOutputEntryCompFuncConst) );
				break;
			}
			//std::cerr << e->Aposition / (this->reader.contigs[e->AcontigID].bases >> 9) << std::endl;
			//assert((e->Aposition / (this->reader.contigs[e->AcontigID].bases >> 9)) < 1024);
			writer << *e;
		}

		/*
		while(true){
			if(!iterators[id]->nextEntry(e))
				break;

			if(!(*e < outQueue.top().data)){
				entry_type hard_copy(e);
				outQueue.push( queue_entry(hard_copy, id, IO::Support::TomahawkOutputEntryCompFuncConst) );
				break;
			}

			writer << *e;
		}
		*/

		/*
		if(iterators[id]->nextEntry(e)){
			// Push new entry into priority queue
			//entry_type hard_copy(e);
			outQueue.push( queue_entry(e, id, IO::Support::TomahawkOutputEntryCompFuncConst) );
		}
		*/


	}

	writer.flush();
	writer.close();

	index_type& index = writer.getIndex();
	std::cerr << index << std::endl;

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
