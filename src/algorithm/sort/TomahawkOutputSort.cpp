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

bool TomahawkOutputSorter::sortMerge(const std::string& inputFile, const std::string& destinationPrefix, const U32 block_size){
	if(!this->reader.Open(inputFile)){
		std::cerr << "failed top open: " << inputFile << std::endl;
		return false;
	}

	if(!this->reader.hasIndex){
		std::cerr << "Does not have TOI index" << std::endl;
		return false;
	}

	this->reader.setWriterType(0);
	this->reader.literals += "\n##tomahawk_sortCommand=" + Helpers::program_string(true);
	if(!this->reader.OpenWriter(destinationPrefix)){
		std::cerr << "failed open" << std::endl;
		return false;
	}

	const U32 n_toi_entries = this->reader.toi_reader.size();
	std::ifstream streams[n_toi_entries];
	IO::TGZFEntryIterator<entry_type>** iterators = new IO::TGZFEntryIterator<entry_type>*[n_toi_entries];

	for(U32 i = 0; i < n_toi_entries; ++i){
		streams[i].open(inputFile);
		streams[i].seekg(this->reader.toi_reader[i].byte_offset);
		iterators[i] = new IO::TGZFEntryIterator<entry_type>(streams[i], 65536, this->reader.toi_reader[i].byte_offset, this->reader.toi_reader[i].byte_offset_end);
	}

	// queue
	queue_type outQueue;

	// draw one from each
	const entry_type* e;
	for(U32 i = 0; i < n_toi_entries; ++i){
		iterators[i]->nextEntry(e);
		entry_type hard_copy(*e); // invoke copy ctor to avoid pointer errors when modifying internal buffer
		outQueue.push( queue_entry(hard_copy, i, IO::Support::TomahawkOutputEntryCompFuncConst) );
	}

	//entry_type prev;
	//prev.AcontigID = 0;
	//prev.Aposition = 0;
	//prev.BcontigID = 0;
	//prev.Bposition = 0;

	// while queue is not empty
	while(outQueue.empty() == false){
		// peek at top entry in queue
		const queue_entry& lowest = outQueue.top();
		const U32 id = lowest.streamID;
		// write the entry from the top of the queue
		//std::cout.write(reinterpret_cast<const char*>(lowest.data), sizeof(entry_type));

		//assert(prev < lowest.data);

		//std::cout << lowest.data << '\n';
		*this->reader.writer << lowest.data;
		// remove this record from the queue
		//std::cerr << prev << '\n' << lowest.data << std::endl << std::endl;
		//prev = lowest.data;

		outQueue.pop();

		// If it is possible to get a new entry from this particular stream
		if(iterators[id]->nextEntry(e)){
			// Push new entry into priority queue
			entry_type hard_copy(*e);
			outQueue.push( queue_entry(hard_copy, id, IO::Support::TomahawkOutputEntryCompFuncConst) );
		}
	}

	// Cleanup
	for(U32 i = 0; i < n_toi_entries; ++i)
		delete iterators[i];
	delete [] iterators;

	return true;
}

}
}
}
