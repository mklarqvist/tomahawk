#ifndef TOTEMPOLE_TOTEMPOLEOUTPUTSORTEDINDEX_H_
#define TOTEMPOLE_TOTEMPOLEOUTPUTSORTEDINDEX_H_

#include <vector>
#include <fstream>

#include "../tomahawk/two/output_entry.h"
#include "../index/index_contig.h"
#include "TotempoleOutputEntry.h"

namespace Tomahawk {
namespace Totempole {

#define TOTEMPOLE_OUTPUT_SORT_CHUNKS 1024
#define TOTEMPOLE_OUTPUT_SORT_SHIFT  9

struct TotempoleOutputSortedIndexBin{
	typedef TotempoleOutputSortedIndexBin self_type;
	typedef TotempoleOutputSortedEntry totempole_entry;

	TotempoleOutputSortedIndexBin(void) :
		n_chunks(0),
		chunks(nullptr)
	{}

	TotempoleOutputSortedIndexBin(const U32 n_chunks) :
		n_chunks(n_chunks),
		chunks(new totempole_entry[n_chunks])
	{}

	~TotempoleOutputSortedIndexBin(){
		delete [] this->chunks;
	}

	bool allocate(const U32 size){
		if(size == 0)
			return false;

		if(this->chunks != nullptr)
			delete [] this->chunks;

		this->n_chunks = size;
		this->chunks = new totempole_entry[size];
		return true;
	}

	totempole_entry& operator[](const U32 p){ return(this->chunks[p]); }
	friend std::ostream& operator<<(std::ostream& stream, const self_type& self){
		for(U32 i = 0; i < self.n_chunks; ++i)
			stream << '\t' << i << '\t' << self.chunks[i] << std::endl;

		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& index){
		stream.write(reinterpret_cast<const char*>(&index.n_chunks), sizeof(U32));
		for(U32 i = 0; i < index.n_chunks; ++i)
			stream << index.chunks[i];

		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& index){
		stream.read(reinterpret_cast<char*>(&index.n_chunks), sizeof(U32));
		index.allocate(index.n_chunks); // allocate memory
		for(U32 i = 0; i < index.n_chunks; ++i)
			stream >> index.chunks[i];

		return(stream);
	}

	U32 n_chunks;
	totempole_entry* chunks;
};

class TotempoleOutputSortedIndex {
	typedef TotempoleOutputSortedIndex self_type;
	typedef TotempoleOutputSortedEntry totempole_entry;
	typedef IO::OutputEntry two_entry;
	typedef Totempole::IndexContigBase contig_type;
	typedef TotempoleOutputSortedIndexBin chunk_type;

public:
	enum TOI_SORTED_ERROR {TOI_SORTED_OK, TOI_SORTED_NO_EXIST, TOI_SORTED_CORRUPTED, TOI_SORTED_INIT};

public:
	TotempoleOutputSortedIndex(const U32 n_contigs, const contig_type* const contigs);
	~TotempoleOutputSortedIndex();

	inline chunk_type& operator[](const U32 p){ return(this->linear_index[p]); }
	inline U32 getBin(const U32& contigID, const U32& position) const{ return(position / ((this->contigs[contigID].n_bases >> TOTEMPOLE_OUTPUT_SORT_SHIFT) + 1)); }
	inline const TOI_SORTED_ERROR& getState(void) const{ return(this->state); }

	void update(const two_entry& entry, const U32& block, const U32& blockOffset);

	bool findOverlap(const U32 contigID, totempole_entry& intervals);
	bool findOverlap(const U32 contigID, const U32 position, totempole_entry& intervals);
	bool findOverlap(const U32 contigID, const U32 from, const U32 to, std::vector<totempole_entry>& intervals);

	friend std::ostream& operator<<(std::ostream& stream, const self_type& self){
		// linear index
		for(U32 i = 0; i < self.n_contigs; ++i)
			stream << "contig: " << i << "\n" << self.linear_index[i] << std::endl;

		// secondary index
		for(U32 i = 0; i < self.n_contigs; ++i)
			stream << self.secondary_index[i] << std::endl;

		return stream;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& index){
		// secondary
		for(U32 i = 0; i < index.n_contigs; ++i)
			stream << index.secondary_index[i];

		// linear
		for(U32 i = 0; i < index.n_contigs; ++i)
			stream << index.linear_index[i];

		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& index){
		// secondary
		for(U32 i = 0; i < index.n_contigs; ++i)
			stream >> index.secondary_index[i];

		// linear
		for(U32 i = 0; i < index.n_contigs; ++i)
			stream >> index.linear_index[i];

		index.state = TOI_SORTED_OK;

		return(stream);
	}

public:
	TOI_SORTED_ERROR state;

private:
	const U32 n_contigs;
	U32 prev_contigIDA;
	U32 prev_chunk;
	const contig_type* const contigs;
	chunk_type* linear_index;
	totempole_entry* secondary_index;
};

} /* namespace Totempole */
} /* namespace Tomahawk */

#endif /* TOTEMPOLE_TOTEMPOLEOUTPUTSORTEDINDEX_H_ */
