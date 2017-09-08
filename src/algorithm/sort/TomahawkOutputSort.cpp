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

	// if index
	// split into thread parts
	// write each


	bool trigger_break = false;
	while(true){
		if(!this->reader.nextBlockUntil(memory_limit)){
			//std::cerr << "failed to get next block: " << this->reader.output_buffer.size() << std::endl;
			trigger_break = true;
		}

		if(this->reader.output_buffer.size() == 0){
			trigger_break = true;
			break;
		}

		assert((this->reader.output_buffer.size() % sizeof(IO::TomahawkOutputEntry)) == 0);

		//std::cerr << Helpers::timestamp("DEBUG") << this->reader.output_buffer.size() << '/' << memory_limit << std::endl;
		std::cerr << Helpers::timestamp("LOG","SORT") << "Sorting: " << Helpers::ToPrettyString(this->reader.output_buffer.size()/sizeof(entry_sort_type)) << " entries" << std::endl;
		std::sort(reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[0]),
				  reinterpret_cast<entry_sort_type*>(&this->reader.output_buffer.data[this->reader.output_buffer.size() - sizeof(entry_sort_type)]));

		//const entry_type* entry;
		//while(this->reader.nextVariantLimited(entry)){
		//	std::cerr << *entry << std::endl;
		//}
		this->reader.writer->write(this->reader.output_buffer);

		std::cerr << Helpers::timestamp("LOG","SORT") << "Writing..." << std::endl;
		if(trigger_break) break;
	}
	std::cerr << this->reader.output_buffer.size() << std::endl;

	this->reader.writer->flush();
	this->reader.writer->close();

	return true;
}

bool TomahawkOutputSorter::sortMerge(const std::string& inputFile){
	if(!this->reader.Open(inputFile)){
		std::cerr << "failed top open: " << inputFile << std::endl;
		return false;
	}

	IO::TGZFControllerStream c;

	//char input_buffer[5012];
	const U32 n_entries_chunk = sizeof(entry_type)*1000;
	BYTE output_buffer[n_entries_chunk];
	//U32 output_buffer_pointer = 0;

	while(this->reader.stream.good()){
		std::cerr << this->reader.stream.tellg() << std::endl;
		if(!c.InflateOpen(this->reader.stream)){
			return false;
		}

		if(!this->reader.stream.good())
			break;

		while(true){
			U32 return_size = 0;
			if(!c.Inflate(this->reader.stream, output_buffer, n_entries_chunk, return_size))
				break;

			//const U32 n_entries = return_size / sizeof(entry_type);
			//const entry_type* const entries = reinterpret_cast<const entry_type* const>(&output_buffer[0]);
			//for(U32 i = 0; i < n_entries; ++i)
			//	std::cout << entries[i] << std::endl;

				//output_buffer_pointer += ret;
			//std::cerr << "done" << std::endl;

			//std::cerr << output_buffer_pointer << '/' << n_entries_chunk << std::endl;
			//std::cerr << "end inner" << std::endl;
			//const U32 remainder = n_entries_chunk % sizeof(entry_type);
			//std::cerr << "remainder: " << remainder << std::endl;
			//memcpy(&output_buffer[0], &output_buffer[n_entries_chunk - remainder], remainder);
			//output_buffer_pointer = remainder;
		}
		this->reader.stream.seekg(IO::Constants::TGZF_BLOCK_FOOTER_LENGTH, std::ios::cur); // tail data
		std::cerr << this->reader.stream.tellg() << '\t' << this->reader.filesize << std::endl;
		if(this->reader.stream.tellg() == this->reader.filesize || !this->reader.stream.good())
			break;
	}

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
