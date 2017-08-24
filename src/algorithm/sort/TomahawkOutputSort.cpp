#include <cassert>

#include "TomahawkOutputSort.h"

namespace Tomahawk{
namespace Algorithm{
namespace Output{

bool TomahawkOutputSorter::sort(const std::string& input, const U64 memory_limit){
	std::vector<std::string> paths = Tomahawk::Helpers::splitLastOf(input, '/', true);
	std::vector<std::string> files = Tomahawk::Helpers::splitLastOf(paths.back(), '.');

	// Todo: if failed to read from file suffix: try to look into file header MAGIC
	if(files[1].size() == 0){
		std::cerr << "could not determine file type from suffix" << std::endl;
		return false;
	}

	std::transform(files[1].begin(), files[1].end(), files[1].begin(), ::tolower);

	// Setup filenames
	std::string outFile = files[0] + "__partial_sort." + Tomahawk::Constants::OUTPUT_LD_SUFFIX;
	std::string outIndex = outFile + '.' + Tomahawk::Constants::OUTPUT_LD_PARTIAL_SORT_INDEX_SUFFIX;
	std::cerr << Helpers::timestamp("LOG","SORT") << "Setting filename to: " << outFile << std::endl;
	std::cerr << Helpers::timestamp("LOG","SORT") << "Setting index filename to: " << outIndex << std::endl;
	outFile = paths[0] + outFile;
	outIndex = paths[0] + outIndex;

	return(this->sort(input, outFile, memory_limit));
}

bool TomahawkOutputSorter::sort(const std::string& input, const std::string& destinationPrefix, const U64 memory_limit){
	const std::string outFile = destinationPrefix + '.' + Tomahawk::Constants::OUTPUT_LD_SUFFIX;
	const std::string outIndex = outFile + '.' + Tomahawk::Constants::OUTPUT_LD_PARTIAL_SORT_INDEX_SUFFIX;

	std::ofstream outIndexStream(outIndex, std::ios::out | std::ios::binary);
	if(!outIndexStream.good()){
		std::cerr << "Faield to open outindex" << std::endl;
		return false;
	}
	outIndexStream.write(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH);
	outIndexStream.write((char*)&Tomahawk::Constants::PROGRAM_VERSION, sizeof(float));

	if(!this->reader.Open(input)){
		std::cerr << "failed top open: " << input << std::endl;
		return false;
	}

	TomahawkOutputSortIndexEntry index;

	//this->reader.writer_output_type = this->reader.WRITER_TYPE::binary;
	this->reader.setWriterType(0);
	if(!this->reader.OpenWriter(outFile)){
		std::cerr << "failed open" << std::endl;
		return false;
	}

	bool trigger_break = false;
	IO::WriterFile& stream = *reinterpret_cast<IO::WriterFile*>(this->reader.writer->getStream());
	std::ofstream& native = stream.getNativeStream();
	index.from = native.tellp();


	while(true){
		if(!this->reader.nextBlockUntil(memory_limit)){
			//std::cerr << "failed to get next block: " << this->reader.output_buffer.size() << std::endl;
			trigger_break = true;
		}

		if(this->reader.output_buffer.size() == 0){
			trigger_break = true;
			break;
		}

		assert((double)this->reader.output_buffer.size()/sizeof(IO::TomahawkOutputEntrySort) == 0);

		std::cerr << Helpers::timestamp("LOG","SORT") << "Sorting: " << Helpers::ToPrettyString(this->reader.output_buffer.size()/sizeof(IO::TomahawkOutputEntrySort)) << " entries" << std::endl;
		std::sort(reinterpret_cast<IO::TomahawkOutputEntrySort*>(&this->reader.output_buffer.data[0]),
				  reinterpret_cast<IO::TomahawkOutputEntrySort*>(&this->reader.output_buffer.data[this->reader.output_buffer.size() - sizeof(IO::TomahawkOutputEntrySort)]));

		//const entry_type* entry;
		//while(this->reader.nextVariantLimited(entry)){
		//	std::cerr << *entry << std::endl;
		//}

		std::cerr << Helpers::timestamp("LOG","SORT") << "Writing..." << std::endl;
		this->reader.writer->write(this->reader.output_buffer);

		index.to = native.tellp();
		std::cerr << "index: " << index << std::endl;
		outIndexStream << index;
		index.from = native.tellp();

		if(trigger_break) break;
	}
	std::cerr << this->reader.output_buffer.size() << std::endl;
	this->reader.writer->writeEOF();
	this->reader.writer->close();

	return true;
}

bool TomahawkOutputSorter::sortMerge(const std::string& inputFile){
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
