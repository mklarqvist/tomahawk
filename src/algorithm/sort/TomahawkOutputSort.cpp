#include <cassert>

#include "TomahawkOutputSort.h"

namespace Tomahawk{
namespace Algorithm{
namespace Output{

bool TomahawkOutputSorter::sort(const std::string& input, const std::string& destinationPrefix, const U64 memory_limit){
	if(!this->reader.Open(input)){
		std::cerr << "failed top open: " << input << std::endl;
		return false;
	}

	this->reader.setWriterType(0);
	this->reader.literals += "\n##tomahawk_sortCommand=" + Helpers::program_string(true);
	this->reader.OpenWriter(destinationPrefix);
	IO::WriterFile toi_writer;
	toi_writer.open(destinationPrefix + ".toi");
	IO::TomahawkOutputSortHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> headIndex(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, this->reader.header.samples, this->reader.header.n_contig);
	toi_writer.getNativeStream() << headIndex;

	// if index
	// split into thread parts
	// write each
	typedef Totempole::TotempoleOutputEntry totempoly_entry;

	// writer
	IO::WriterFile& stream = *reinterpret_cast<IO::WriterFile*>(this->reader.writer->getStream());

	bool trigger_break = false;
	while(true){
		std::cerr << Helpers::timestamp("LOG","SORT") << "Reading..." << std::endl;
		if(!this->reader.nextBlockUntil(memory_limit)){
			//std::cerr << "failed to get next block: " << this->reader.output_buffer.size() << std::endl;
			trigger_break = true;
		}

		if(this->reader.output_buffer.size() == 0){
			trigger_break = true;
			break;
		}

		assert((this->reader.output_buffer.size() % sizeof(entry_type)) == 0);

		//std::cerr << Helpers::timestamp("DEBUG") << this->reader.output_buffer.size() << '/' << memory_limit << std::endl;
		std::cerr << Helpers::timestamp("LOG","SORT") << "Sorting: " << Helpers::ToPrettyString(this->reader.output_buffer.size()/sizeof(entry_sort_type)) << " entries" << std::endl;
		std::sort(reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[0]),
				  reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[this->reader.output_buffer.size()]));


		std::cerr << Helpers::timestamp("LOG","SORT") << "Indexing..." << std::endl;
		const entry_type* entry;
		if(!this->reader.nextVariantLimited(entry)){
			std::cerr << "no data" << std::endl;
			return false;
		}
		totempoly_entry totempole;
		totempole.contigIDA = entry->AcontigID;
		totempole.contigIDB = entry->BcontigID;
		totempole.minPositionA = entry->Aposition;
		totempole.minPositionB = entry->Bposition;
		totempole.entries = 1;
		totempole.byte_offset = stream.getNativeStream().tellp();
		totempole.uncompressed_size = this->reader.output_buffer.size();

		const entry_type* prev;
		std::swap(prev, entry);

		while(this->reader.nextVariantLimited(entry)){
			//std::cerr << *entry << std::endl;
			if(totempole.contigIDA != -1){
				if(entry->Aposition < prev->Aposition || entry->Aposition < totempole.minPositionA){
					totempole.minPositionA = -1;
					totempole.maxPositionA = -1;
				} else if(totempole.minPositionA != -1) totempole.maxPositionA = entry->Aposition;

				if(totempole.contigIDA != entry->AcontigID){
					totempole.contigIDA = -1;
					totempole.minPositionA = -1;
					totempole.maxPositionA = -1;
				}
			}

			if(totempole.contigIDB != -1){
				if(entry->Bposition < prev->Bposition || entry->Bposition < totempole.minPositionB){
					totempole.minPositionB = -1;
					totempole.maxPositionB = -1;
				} else if(totempole.minPositionB != -1) totempole.maxPositionB = entry->Bposition;

				if(totempole.contigIDB != entry->BcontigID){
					totempole.contigIDB = -1;
					totempole.minPositionB = -1;
					totempole.maxPositionB = -1;
				}
			}
			++totempole.entries;
			std::swap(prev, entry);
		}

		this->reader.writer->write(this->reader.output_buffer);
		totempole.byte_offset_end = stream.getNativeStream().tellp();
		std::cerr << totempole << std::endl;
		toi_writer.getNativeStream() << totempole;

		std::cerr << Helpers::timestamp("LOG","SORT") << "Writing..." << std::endl;
		if(trigger_break) break;
	}
	toi_writer.flush();
	toi_writer.close();

	this->reader.writer->flush();
	this->reader.writer->close();

	return true;
}

bool TomahawkOutputSorter::sortMerge(const std::string& inputFile){
	if(!this->reader.Open(inputFile)){
		std::cerr << "failed top open: " << inputFile << std::endl;
		return false;
	}

	if(!this->reader.hasIndex){
		std::cerr << "Does not have TOI index" << std::endl;
		return false;
	}

	const U32 n_toi_entries = this->reader.toi_reader.size();
	std::ifstream streams[n_toi_entries];
	IO::TGZFEntryIterator<entry_type>** iterators = new IO::TGZFEntryIterator<entry_type>*[n_toi_entries];

	for(U32 i = 0; i < n_toi_entries; ++i){
		std::cerr << i << '\t' << this->reader.toi_reader[i] << std::endl;
		streams[i].open(inputFile);
		streams[i].seekg(this->reader.toi_reader[i].byte_offset);
		iterators[i] = new IO::TGZFEntryIterator<entry_type>(streams[i], 65536, this->reader.toi_reader[i].byte_offset, this->reader.toi_reader[i].byte_offset_end);
	}

	//IO::TGZFEntryIterator<entry_type> c(this->reader.stream, 65536, 0, this->reader.filesize);
	/*
	const entry_type* entry;

	U64 count = 0;

	for(U32 i = 0; i < n_toi_entries; ++i){
		while(iterators[i]->nextEntry(entry)){
			++count;
		}
		std::cerr << "done " << i << std::endl;
	}

	//while(c.nextEntry(entry)){
	//	++count;
	//}
	std::cerr << "count: " << count << std::endl;



	return true;
	*/

	// queue
	queue_type outQueue;

	// draw one from each
	const entry_type* e;
	for(U32 i = 0; i < n_toi_entries; ++i){
		iterators[i]->nextEntry(e);
		entry_type hard_copy(*e);
		outQueue.push( queue_entry(hard_copy, i, IO::Support::TomahawkOutputEntryCompFuncConst) );
	}

	entry_type prev;
	prev.AcontigID = 0;
	prev.Aposition = 0;
	prev.BcontigID = 0;
	prev.Bposition = 0;

	// while queue is not empty
	while(outQueue.empty() == false){
		// peek at top entry in queue
		const queue_entry& lowest = outQueue.top();
		const U32 id = lowest.streamID;
		// write the entry from the top of the queue
		//std::cout.write(reinterpret_cast<const char*>(lowest.data), sizeof(entry_type));

		assert(prev < lowest.data);

		//std::cout << lowest.data << '\n';
		// remove this record from the queue
		std::cerr << prev << '\n' << lowest.data << std::endl << std::endl;
		prev = lowest.data;

		outQueue.pop();

		// If it is possible to get a new entry from this particular stream
		if(iterators[id]->nextEntry(e)){
			// Push new entry into priority queue
			entry_type hard_copy(*e);
			outQueue.push( queue_entry(hard_copy, id, IO::Support::TomahawkOutputEntryCompFuncConst) );
		}
	}

	for(U32 i = 0; i < n_toi_entries; ++i)
		delete iterators[i];
	delete [] iterators;

	return true;

	/*
	std::cerr << "attempting to open: " << inputFile + '.' + Tomahawk::Constants::OUTPUT_LD_PARTIAL_SORT_INDEX_SUFFIX << std::endl;
	std::ifstream indexStream(inputFile + '.' + Tomahawk::Constants::OUTPUT_LD_PARTIAL_SORT_INDEX_SUFFIX, std::ios::binary | std::ios::in | std::ios::ate);
	if(!indexStream.good()){
		std::cerr << "bad index stream" << std::endl;
		return false;
	}
	const U64 filesize_index = indexStream.tellg();
	indexStream.seekg(0);

	// parse header
	partial_header_type head;
	indexStream >> head;
	if(!head.validate()){
		std::cerr << "Failed to validate header" << std::endl;
		return false;
	}
	std::cerr << std::string(&head.header[0], Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH) << '\t' << head.version << std::endl;
	std::cerr << "position now: " << indexStream.tellg() << '/' << filesize_index << std::endl;
	std::cerr << "entries: " << (filesize_index - indexStream.tellg())/(2*sizeof(U64)) << std::endl;

	char* indexHeaderEntries = new char[filesize_index - indexStream.tellg()];
	const U32 indexEntryEnd = (filesize_index - indexStream.tellg())/(2*sizeof(U64));
	indexStream.read(indexHeaderEntries, filesize_index - indexStream.tellg());
	const partial_header_entry_type* const entries = reinterpret_cast<const partial_header_entry_type* const>(&indexHeaderEntries[0]);

	sort_reader* sortEntries = new sort_reader[indexEntryEnd];

	std::cerr << Helpers::timestamp("LOG", "SORT") << "Opening " << indexEntryEnd << " file handles..." << std::endl;
	for(U32 i = 0; i < indexEntryEnd; ++i){
		//std::cerr << entries[i] << std::endl;

		if(!sortEntries[i].setup(inputFile, entries[i].from, entries[i].to, 1000000 - (1000000 % sizeof(entry_type)))){
			std::cerr << "failed setup" << std::endl;
		}
	}

	//
	if(!this->kwayMerge(sortEntries, indexEntryEnd, std::cout)){
		delete [] sortEntries;
		delete [] indexHeaderEntries;
		return false;
	}

	delete [] sortEntries;
	delete [] indexHeaderEntries;

	return true;
	*/
}

}
}
}
