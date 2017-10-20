#ifndef TOTEMPOLE_TOTEMPOLEOUTPUTINDEXWRITER_H_
#define TOTEMPOLE_TOTEMPOLEOUTPUTINDEXWRITER_H_

namespace Tomahawk{
namespace Totempole{

class TomahawkOutputWriterIndex{
	typedef TomahawkOutputWriterIndex self_type;
	typedef IO::TomahawkOutputEntry entry_type;
	typedef TotempoleOutputEntry totempole_output_entry;
	typedef Tomahawk::IO::TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC_LENGTH> toi_header_type;
	typedef Tomahawk::IO::TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> two_header_type;
	typedef Totempole::TotempoleOutputSortedIndex index_type;
	typedef Totempole::TotempoleContigBase contig_type;

public:
	TomahawkOutputWriterIndex(const contig_type* contigs, const U32& n_contigs, const two_header_type& header):
		currentBlockID(0),
		currentOffset(0),
		index(n_contigs, contigs),
		header(Tomahawk::Constants::WRITE_HEADER_LD_SORT_MAGIC, header.samples, header.n_contig)
	{
		this->header.controller = header.controller;
		this->header.n_entries = header.n_entries;
	}
	~TomahawkOutputWriterIndex(){}

	bool Open(const std::string filename){
		this->stream.open(filename, std::ios::out | std::ios::binary);
		if(!this->stream.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TOI") << "Failed to open " << filename << "..." << std::endl;
			return false;
		}
		std::cerr << Tomahawk::Helpers::timestamp("LOG", "TOI") << "Opening: " << filename << "..." << std::endl;
		return true;
	}

	bool WriteHeader(void){
		this->stream << this->header;

		return true;
	}

	// Not expanded and not sorted
	// Basic TWI entries
	void flushBlock(const U64& fromOffset, const U64& toOffset, const U32& uncompressed_size, const U32& n_entries){
		this->totempole_entry.entries = n_entries;
		this->totempole_entry.byte_offset = fromOffset;
		this->totempole_entry.byte_offset_end = toOffset;
		this->totempole_entry.uncompressed_size = uncompressed_size;
		this->stream << this->totempole_entry;
		++this->currentBlockID;
	}

	// If the TWO file is expanded and sorted
	// then make use of the index
	bool UpdateIndexed(const entry_type& entry){
		++this->totempole_entry.entries;
		this->index.update(entry, this->currentBlockID, this->currentOffset);
		++this->currentOffset;

		return true;
	}

	void flushBlockIndexed(const U64& fromOffset, const U64& toOffset, const U32& uncompressed_size){
		this->totempole_entry.byte_offset = fromOffset;
		this->totempole_entry.byte_offset_end = toOffset;
		this->totempole_entry.uncompressed_size = uncompressed_size;
		this->stream << this->totempole_entry;
		this->totempole_entry.reset();
		++this->currentBlockID;
		this->currentOffset = 0;
	}

	inline void writeIndex(void){ this->stream << this->index; }

	// Set and get
	inline toi_header_type& getHeader(void){ return(this->header); }
	inline index_type& getIndex(void){ return(this->index); }
	inline const U32& blocksWritten(void) const{ return(this->currentBlockID); }

private:
	U32 currentBlockID;
	U32 currentOffset;
	totempole_output_entry totempole_entry;
	index_type index;
	std::ofstream stream;
	toi_header_type header;
};

}
}

#endif /* TOTEMPOLE_TOTEMPOLEOUTPUTINDEXWRITER_H_ */
