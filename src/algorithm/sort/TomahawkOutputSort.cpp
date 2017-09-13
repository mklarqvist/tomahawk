#include <cassert>

#include "TomahawkOutputSort.h"
#include "../../io/compression/TGZFEntryIterator.h"

namespace Tomahawk{
namespace Algorithm{
namespace Output{

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

	//
	this->reader.setWriterType(0);
	this->reader.literals += "\n##tomahawk_partialSortCommand=" + Helpers::program_string(true);
	this->reader.OpenWriter(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX);
	IO::WriterFile toi_writer;
	toi_writer.open(basePath + baseName + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX + '.' + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX);
	IO::TomahawkOutputSortHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> headIndex(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, this->reader.header.samples, this->reader.header.n_contig);
	toi_writer.getNativeStream() << headIndex;

	// if index
	// split into thread parts
	// write each

	// writer
	IO::WriterFile& stream = *reinterpret_cast<IO::WriterFile*>(this->reader.writer->getStream());

	bool trigger_break = false;
	while(true){
		std::cerr << Helpers::timestamp("LOG","SORT") << "Reading..." << std::endl;
		if(!this->reader.nextBlockUntil(memory_limit))
			trigger_break = true;

		if(this->reader.output_buffer.size() == 0){
			trigger_break = true;
			break;
		}

		assert((this->reader.output_buffer.size() % sizeof(entry_type)) == 0);

		std::cerr << Helpers::timestamp("LOG","SORT") << "Sorting: " << Helpers::ToPrettyString(this->reader.output_buffer.size()/sizeof(entry_sort_type)) << " entries" << std::endl;
		std::sort(reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[0]),
				  reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[this->reader.output_buffer.size()]));


		std::cerr << Helpers::timestamp("LOG","SORT") << "Indexing..." << std::endl;
		const entry_type* entry;
		if(!this->reader.nextVariantLimited(entry)){
			std::cerr << Helpers::timestamp("ERROR","SORT") << "No data" << std::endl;
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
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed to open: " << inputFile << "..." << std::endl;
		return false;
	}

	if(!this->reader.hasIndex){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Does not have TOI index" << std::endl;
		return false;
	}

	this->reader.literals += "\n##tomahawk_mergeSortCommand=" + Helpers::program_string(true);

	writer_type writer(this->reader.contigs, &this->reader.header);
	if(!writer.open(destinationPrefix)){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "Failed open" << std::endl;
		return false;
	}
	writer.writeHeader(this->reader.literals);

	const U32 n_toi_entries = this->reader.toi_reader.size();
	std::ifstream* streams = new std::ifstream[n_toi_entries];
	IO::TGZFEntryIterator<entry_type>** iterators = new IO::TGZFEntryIterator<entry_type>*[n_toi_entries];

	std::cerr << Helpers::timestamp("LOG", "SORT") << "Opening " << n_toi_entries << " file handles...";
	for(U32 i = 0; i < n_toi_entries; ++i){
		streams[i].open(inputFile);
		streams[i].seekg(this->reader.toi_reader[i].byte_offset);
		iterators[i] = new IO::TGZFEntryIterator<entry_type>(streams[i], 65536, this->reader.toi_reader[i].byte_offset, this->reader.toi_reader[i].byte_offset_end);
	}
	std::cerr << " Done!" << std::endl;

	// queue
	queue_type outQueue;

	//
	std::cerr << Helpers::timestamp("LOG", "SORT") << "Merging..." << std::endl;

	// draw one from each
	const entry_type* e;
	for(U32 i = 0; i < n_toi_entries; ++i){
		iterators[i]->nextEntry(e);
		entry_type hard_copy(e); // invoke copy ctor to avoid pointer errors when modifying internal buffer
		outQueue.push( queue_entry(hard_copy, i, IO::Support::TomahawkOutputEntryCompFuncConst) );
	}

	if(outQueue.empty()){
		std::cerr << Helpers::timestamp("ERROR","SORT") << "No data in queue..." << std::endl;
		return false;
	}

	// index
	writer.setPrevEntry(outQueue.top().data);
	writer.setPrevEntryFirst(outQueue.top().data);

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
		if(iterators[id]->nextEntry(e)){
			// Push new entry into priority queue
			entry_type hard_copy(e);
			outQueue.push( queue_entry(hard_copy, id, IO::Support::TomahawkOutputEntryCompFuncConst) );
		}
	}
	//std::cerr << "Total output entries: " << outputEntries << std::endl;

	writer.flush();
	writer.close();

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
