#include "TotempoleOutputSortedIndex.h"

namespace Tomahawk {
namespace Totempole {

TotempoleOutputSortedIndex::TotempoleOutputSortedIndex(const U32 n_contigs, const contig_type* const contigs) :
	state(TOI_SORTED_INIT),
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

TotempoleOutputSortedIndex::~TotempoleOutputSortedIndex(){
	delete [] this->linear_index;
	delete [] this->secondary_index;
}


void TotempoleOutputSortedIndex::update(const two_entry& entry, const U32& block, const U32& blockOffset){
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
	}

	// secondary index
	// in case there is no switch then last entry is +1 (non-closed interval)
	this->linear_index[entry.AcontigID][chunk].update(block, blockOffset + 1);
	this->secondary_index[entry.AcontigID].update(block, blockOffset + 1);
	this->prev_contigIDA = entry.AcontigID;
	this->prev_chunk = chunk;
}


bool TotempoleOutputSortedIndex::findOverlap(const U32 contigID, totempole_entry& intervals){
	if(this->secondary_index[contigID].fromBlock == -1)
		return false;

	intervals = this->secondary_index[contigID];

	return true;
}

bool TotempoleOutputSortedIndex::findOverlap(const U32 contigID, const U32 position, totempole_entry& intervals){
	if(this->secondary_index[contigID].fromBlock == -1)
		return false;

	const U32 chunk = position / (this->contigs[contigID].bases >> TOTEMPOLE_OUTPUT_SORT_SHIFT);

	// Impossible chunk
	if(chunk >= TOTEMPOLE_OUTPUT_SORT_CHUNKS)
		return false;

	if(this->linear_index[contigID][chunk].fromBlock == -1)
		return false;

	intervals = this->linear_index[contigID][chunk];
	return true;
}

bool TotempoleOutputSortedIndex::findOverlap(const U32 contigID, const U32 from, const U32 to, std::vector<totempole_entry>& intervals){
	if(this->secondary_index[contigID].fromBlock == -1)
		return false;

	const U32 chunkFrom = from / (this->contigs[contigID].bases >> TOTEMPOLE_OUTPUT_SORT_SHIFT);
	U32 chunkTo         = to   / (this->contigs[contigID].bases >> TOTEMPOLE_OUTPUT_SORT_SHIFT);

	if(chunkFrom >= TOTEMPOLE_OUTPUT_SORT_CHUNKS)
		return false;

	if(chunkTo >= TOTEMPOLE_OUTPUT_SORT_CHUNKS)
		chunkTo = TOTEMPOLE_OUTPUT_SORT_CHUNKS - 1;

	// Special case when both values fall into the same chunk
	if(chunkFrom == chunkTo){
		if(this->linear_index[contigID][chunkFrom].fromBlock == -1)
			return false;

		intervals.push_back(this->linear_index[contigID][chunkFrom]);
		return true;
	}

	totempole_entry entry;
	U32 i = chunkFrom;
	U32 prevHit = 0;
	// Get first
	for(; i < chunkTo; ++i){
		if(this->linear_index[contigID][i].fromBlock != -1){
			entry = this->linear_index[contigID][i];
			prevHit = i;
			break;
		}
	}
	++i;

	// Continue until end
	for(; i < chunkTo; ++i){
		if(this->linear_index[contigID][i].fromBlock != -1){
			// extend
			if(i - prevHit == 1){
				entry.toBlock = this->linear_index[contigID][i].toBlock;
				entry.toBlock_entries_offset = this->linear_index[contigID][i].toBlock_entries_offset;
				prevHit = i;
			} else {
				intervals.push_back(entry);
				entry = this->linear_index[contigID][i];
			}
		}
	}

	// push back last
	if(i - prevHit == 1){
		entry.toBlock = this->linear_index[contigID][i].toBlock;
		entry.toBlock_entries_offset = this->linear_index[contigID][i].toBlock_entries_offset;
		intervals.push_back(entry);
	}

	return(intervals.size() > 0);
}


} /* namespace Totempole */
} /* namespace Tomahawk */
