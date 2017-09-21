#ifndef TOTEMPOLE_TOTEMPOLEOUTPUTSORTEDINDEX_H_
#define TOTEMPOLE_TOTEMPOLEOUTPUTSORTEDINDEX_H_

#include <vector>
#include <fstream>

#include "../tomahawk/TomahawkOutput/TomahawkOutputEntry.h"
#include "TotempoleContig.h"
#include "TotempoleOutputEntry.h"

namespace Tomahawk {
namespace Totempole {

#define TOTEMPOLE_OUTPUT_SORT_CHUNKS 1024
#define TOTEMPOLE_OUTPUT_SORT_SHIFT 9

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

	U32 n_chunks;
	totempole_entry* chunks;
};

class TotempoleOutputSortedIndex {
	typedef TotempoleOutputSortedIndex self_type;
	typedef TotempoleOutputSortedEntry totempole_entry;
	typedef IO::TomahawkOutputEntry two_entry;
	typedef Totempole::TotempoleContigBase contig_type;
	typedef TotempoleOutputSortedIndexBin chunk_type;

public:
	TotempoleOutputSortedIndex(const U32 n_contigs, const contig_type* const contigs) :
		n_contigs(n_contigs),
		prev_contigIDA(n_contigs + 1), // init to impossible id
		prev_chunk(0),
		contigs(contigs),
		linear_index(new chunk_type[n_contigs]),
		secondary_index(new totempole_entry[n_contigs])
	{
		for(U32 i = 0; i < this->n_contigs; ++i){
			this->linear_index[i].allocate(TOTEMPOLE_OUTPUT_SORT_CHUNKS);
		}
	}

	~TotempoleOutputSortedIndex(){
		delete [] this->linear_index;
		delete [] this->secondary_index;
	}

	bool parse(std::ifstream& stream);
	chunk_type& operator[](const U32 p){ return(this->linear_index[p]); }

	inline void update(const two_entry& entry, const U32& block, const U32& blockOffset){
		const U32 chunk = this->getBin(entry.AcontigID, entry.Aposition);

		// Switch in chromosome
		if(entry.AcontigID != this->prev_contigIDA){
			// This is only ever NOT TRUE for the first entry
			if(prev_contigIDA != this->n_contigs + 1){
				// Update previous block data with end positions
				totempole_entry& temp = this->linear_index[this->prev_contigIDA][this->prev_chunk];
				temp.toBlock = block;
				temp.toBlock_entries_offset = blockOffset;

				// Update current block data with start positions
				totempole_entry& temp2 = this->linear_index[entry.AcontigID][chunk];
				temp2.fromBlock = block;
				temp2.fromBlock_entries_offset = blockOffset;

				// Secondary
				this->secondary_index[this->prev_contigIDA].toBlock = block;
				this->secondary_index[this->prev_contigIDA].toBlock_entries_offset = blockOffset;
				this->secondary_index[entry.AcontigID].fromBlock = block;
				this->secondary_index[entry.AcontigID].fromBlock_entries_offset = blockOffset;
				std::cerr << "switch in contigID: " << this->prev_contigIDA << "->" << entry.AcontigID << std::endl;
			} else {
				// For the first entry
				const U32 chunk = this->getBin(entry.AcontigID, entry.Aposition);
				totempole_entry& temp = this->linear_index[entry.AcontigID][chunk];
				temp.fromBlock = block;
				temp.fromBlock_entries_offset = blockOffset;

				this->secondary_index[entry.AcontigID].fromBlock = block;
				this->secondary_index[entry.AcontigID].fromBlock_entries_offset = blockOffset;
			}
		}
		// Switch in chunk
		else if(chunk != this->prev_chunk){
			this->linear_index[this->prev_contigIDA][this->prev_chunk].toBlock = block;
			this->linear_index[this->prev_contigIDA][this->prev_chunk].toBlock_entries_offset = blockOffset;

			this->linear_index[entry.AcontigID][chunk].fromBlock = block;
			this->linear_index[entry.AcontigID][chunk].fromBlock_entries_offset = blockOffset;
			std::cerr << "switch in chunk: " << this->prev_chunk << "->" << chunk << std::endl;
		}

		// secondary index
		// in case there is no switch then last entry is +1 (non-closed interval)
		this->linear_index[entry.AcontigID][chunk].update(block, blockOffset + 1);
		this->secondary_index[entry.AcontigID].update(block, blockOffset + 1);
		this->prev_contigIDA = entry.AcontigID;
		//this->prev_contigIDB = entry.BcontigID;
		this->prev_chunk = chunk;
	}

	inline U32 getBin(const U32& contigID, const U32& position) const{ return(position / (this->contigs[contigID].bases >> TOTEMPOLE_OUTPUT_SORT_SHIFT)); }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& self){
		for(U32 i = 0; i < 2; ++i)
			stream << "contig: " << i << "\n" << self.linear_index[i] << std::endl;

		for(U32 i = 0; i < self.n_contigs; ++i)
			stream << self.secondary_index[i] << std::endl;

		return stream;
	}

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
