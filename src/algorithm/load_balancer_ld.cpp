#include "load_balancer_ld.h"

namespace tomahawk{

LoadBalancerLD::LoadBalancerLD() :
		selected_chunk(0),
		n_desired_chunks(1),
		n_comparisons_chunk(0)
{}

LoadBalancerLD::~LoadBalancerLD(){}

bool LoadBalancerLD::getSelectedLoad(){
	const value_type& selected = this->blocks[this->selected_chunk];

	// If there are both equal
	if(selected.fromRow == selected.fromColumn && selected.toRow == selected.toColumn){
		this->data_to_load.push_back(std::pair<U32, U32>(selected.fromRow, selected.toRow));
	} else {
		// No overlap
		this->data_to_load.push_back(std::pair<U32, U32>(selected.fromRow,    selected.toRow));
		this->data_to_load.push_back(std::pair<U32, U32>(selected.fromColumn, selected.toColumn));
	}

	return true;
}

bool LoadBalancerLD::getSelectedLoadInterval(const std::vector< std::pair<U32,U32> >& intervals){
	if(intervals.size() == 0) return false;
	// Transform absolute offsets into relative offsets
	// and store the absolute block offsets
	// e.g.
	// 1: 0->6
	// 2: 265->450
	// 3: 530->559
	// Cumulative distribution = [0,6,191,220]
	// Resulting relative offsets
	// 1: 0->6
	// 2: 6->191
	// 3: 191->220
	//
	// Therefore: to find the overlap of [A,B] and [X,Y] involves
	// simple interval intersections.
	std::vector< std::pair<U32, U32> > old = this->data_to_load; // copy old
	this->data_to_load.clear(); // clear old

	std::vector<U32> cumsum(1, 0); // start with 0
	std::cerr << cumsum[0] << std::endl;
	for(U32 i = 0; i < intervals.size(); ++i){
		cumsum.push_back(intervals[i].second - intervals[i].first + cumsum.back());
		std::cerr << cumsum.back() << std::endl;
	}

	std::vector< std::pair<U32, U32> > relative_intervals;
	for(U32 i = 0; i < intervals.size(); ++i){
		relative_intervals.push_back(std::pair<U32,U32>(cumsum[i], cumsum[i+1]));
		std::cerr << "Relative: " << relative_intervals.back().first << "->" << relative_intervals.back().second << std::endl;
		std::cerr << "Absolute: " << intervals[i].first << "->" << intervals[i].second << std::endl;
	}

	std::cerr << "Start desired" << std::endl;
	for(U32 i = 0 ; i < old.size(); ++i){
		std::cerr << old[i].first << "\t" << old[i].second << std::endl;
	}
	std::cerr << this->blocks[this->selected_chunk] << std::endl;

	std::cerr << "Start map" << std::endl;

	// Cycle over desired relative regions
	// overlap relative vs absolute
	// [a, b] overlaps with [x, y] iff b > x and a < y.
	// a = old[i].first --- desired relative
	// b = old[i].second --- desired relative
	// x = relative_intervals[j].first --- relative
	// y = relative_intervals[j].second --- relative
	for(U32 i = 0; i < old.size(); ++i){ // cycle over relative data to load
		for(U32 j = 0; j < relative_intervals.size(); ++j){ // cycle over relative intervals
			if(old[i].second > relative_intervals[j].first && old[i].first < relative_intervals[j].second){
				std::cerr << "desired collapsed: " << old[i].first << "-" << old[i].second << " overlapping with relative interval " << relative_intervals[j].first << "-" << relative_intervals[j].second << std::endl;
				std::cerr << "absolute: " << intervals[j].first << "->" << intervals[j].second << std::endl;

				U32 start = 0;
				U32 end   = 0;

				// if A < X -> start = X
				if(old[i].first < relative_intervals[j].first) start = relative_intervals[j].first;
				// else start = A
				else start = old[i].first;
				// if B > Y -> end = Y
				if(old[i].second > relative_intervals[j].second) end = relative_intervals[j].second;
				// else end = B
				else end = old[i].second;

				std::cerr << start << "->" << end << std::endl;

				// Remove start position from these
				std::cerr << start - relative_intervals[j].first << "->" << end - relative_intervals[j].first << std::endl;
				std::cerr << "final: " << intervals[j].first + start - relative_intervals[j].first << "->" << intervals[j].first + end - relative_intervals[j].first << std::endl;
				this->data_to_load.push_back(std::pair<U32,U32>(intervals[j].first + start - relative_intervals[j].first, intervals[j].first + end - relative_intervals[j].first));

			} else {

			}
		}
	}

	//exit(1);

	// Dedupe
	std::sort(this->data_to_load.begin(), this->data_to_load.end());
	std::vector<std::pair<U32,U32>> deduped;
	deduped.push_back(this->data_to_load[0]);
	for(U32 i = 1; i < this->data_to_load.size(); ++i){
		if(this->data_to_load[i].first <= deduped.back().second){
			if(deduped.back().second < this->data_to_load[i].second)
				deduped.back().second = this->data_to_load[i].second;
		}
		else deduped.push_back(this->data_to_load[i]);
	}

	for(U32 i = 0; i < this->data_to_load.size(); ++i){
		std::cerr << "original: " << this->data_to_load[i].first << "->" << this->data_to_load[i].second << std::endl;
	}

	for(U32 i = 0; i < deduped.size(); ++i){
		std::cerr << "deduped: " << deduped[i].first << "->" << deduped[i].second << std::endl;
	}
	this->data_to_load = deduped;
	return true;
}

bool LoadBalancerLD::getSelectedLoadThreadsInterval(const reader_type& reader, const U32 threads, const std::vector< std::pair<U32,U32> >& intervals){
	std::cerr << this->selected_chunk << "/" << this->blocks.size() << std::endl;
	const value_type& selected = this->blocks[this->selected_chunk];
	std::cerr << selected << std::endl;
	this->thread_distribution.resize(threads);

	// If a square is attached to the diagonal
	// In this case, workload has to be distrubuted on the
	// upper triangular of comparisons to avoid duplicate
	// computation
	if(selected.isDiagonal()){
		std::cerr << "diagonal" << std::endl;
		U64 n_comparisons_thread = selected.getSize() / threads;
		if(threads == 1) n_comparisons_thread = std::numeric_limits<U64>::max();

		U64 n_blocks_incrementor = 0;
		std::pair<U32, U32> current;
		U32 currentThread = 0;

		for(U32 i = selected.fromRow; i < selected.toRow; ++i){
			current.first  = i;
			current.second = i;

			for(U32 j = i; j < selected.toColumn; ++j){
				++n_blocks_incrementor;

				if(n_blocks_incrementor > n_comparisons_thread){
					current.second = j + 1; // not inclusive range [A, B)
					if(current.second > selected.toColumn) current.second = selected.toColumn;
					//std::cerr << "break counter: " << current.first << "->" << current.second << " with " << n_variants_counter << "/" << n_variants_partition << " of " << this->n_comparisons_chunk << std::endl;
					this->thread_distribution[currentThread].push_back(LoadBalancerThread(i - selected.fromRow, current.first - selected.fromRow, current.second - selected.fromRow));

					current.first  = j + 1;
					current.second = j;
					n_blocks_incrementor = 0;
					++currentThread;
				}
			}

			//std::cerr << current.second << std::endl;
			if(current.second != selected.toRow){
				current.second = selected.toRow;
				//std::cerr << "break new line: " << current.first << "->" << current.second << " with " << n_variants_counter <<"/" << n_variants_partition << std::endl;
				this->thread_distribution[currentThread].push_back(LoadBalancerThread(i - selected.fromRow, current.first - selected.fromRow, current.second - selected.fromRow));
				//parts.push_back(current);
				if(n_blocks_incrementor > n_comparisons_thread){
					//std::cerr << "break here" << std::endl;
					n_blocks_incrementor = 0;
					++currentThread;
				}
			}
		}

		std::cerr << "currentThread: " << currentThread << std::endl;
		std::cerr << selected.toRow - selected.fromRow << std::endl;
		for(U32 i = 0; i < this->thread_distribution.size(); ++i){
			for(U32 j = 0; j < this->thread_distribution[i].size(); ++j){
				std::cerr << "idx: "<< i << "; row: " << this->thread_distribution[i][j].row << ":" << this->thread_distribution[i][j].fromColumn << "->" << this->thread_distribution[i][j].toColumn << "/" << reader.getIndex().getContainer().size() << std::endl;
			}
		}
		std::cerr << "here" << std::endl;


	} else {
		std::cerr << "square" << std::endl;
		const value_type& selectedRow = this->blocks[this->selected_chunk];
		const value_type& selectedCol = this->blocks[this->selected_chunk];

		U64 n_comparisons_thread = selected.getSize() / threads;
		if(threads == 1) n_comparisons_thread = std::numeric_limits<U64>::max();

		U64 n_blocks_incrementor = 0;
		std::pair<U32, U32> current;
		U32 currentThread = 0;


		for(U32 i = selectedRow.fromRow; i < selectedRow.toRow; ++i){
			current.first  = selectedCol.fromColumn;
			current.second = selectedCol.fromColumn;

			for(U32 j = selectedCol.fromColumn; j < selectedCol.toColumn; ++j){
				++n_blocks_incrementor;

				if(n_blocks_incrementor > n_comparisons_thread){
					current.second = j + 1; // not n_comparisons_thread range [A, B)
					if(current.second > selected.toColumn) current.second = selected.toColumn;
					//std::cerr << "break counter: " << current.first << "->" << current.second << " with " << n_variants_counter << "/" << n_variants_partition << std::endl;
					this->thread_distribution[currentThread].push_back(LoadBalancerThread(
							i - selectedRow.fromRow,
							(selectedRow.toRow - selectedRow.fromRow) + current.first - selectedCol.fromColumn,
							(selectedRow.toRow - selectedRow.fromRow) + current.second - selectedCol.fromColumn));

					current.first = j + 1;
					current.second = j;
					n_blocks_incrementor = 0;
					++currentThread;
				}
			}

			//std::cerr << current.second << std::endl;
			if(current.second != selected.toColumn){
				current.second = selected.toColumn;
				//std::cerr << "break new line: " << current.first << "->" << current.second << std::endl;
				this->thread_distribution[currentThread].push_back(LoadBalancerThread(
						i - selectedRow.fromRow,
						(selectedRow.toRow - selectedRow.fromRow) + current.first - selectedCol.fromColumn,
						(selectedRow.toRow - selectedRow.fromRow) + current.second - selectedCol.fromColumn));
				//parts.push_back(current);
				if(n_blocks_incrementor > n_comparisons_thread){
					n_blocks_incrementor = 0;
					++currentThread;
				}
			}
		}

		std::cerr << "currentThread: " << currentThread << std::endl;
		std::cerr << selected.toRow - selected.fromRow << std::endl;
		for(U32 i = 0; i < this->thread_distribution.size(); ++i){
			for(U32 j = 0; j < this->thread_distribution[i].size(); ++j){
				std::cerr << "idx: "<< i << "; row: " << this->thread_distribution[i][j].row << ":" << this->thread_distribution[i][j].fromColumn << "->" << this->thread_distribution[i][j].toColumn << "/" << reader.getIndex().getContainer().size() << std::endl;
			}
		}
		std::cerr << "here" << std::endl;

	}

	//for(U32 i = 0; i < this->data_to_load.size(); ++i)
	//	std::cerr << i << '\t' << this->data_to_load[i].first << "->" << this->data_to_load[i].second << std::endl;

	std::cerr << "returning" << std::endl;
	return true;
}

bool LoadBalancerLD::getSelectedLoadThreads(const reader_type& reader, const U32 threads){
	const value_type& selected = this->blocks[this->selected_chunk];
	this->thread_distribution.resize(threads);

	// If a square is attached to the diagonal
	// In this case, workload has to be distrubuted on the
	// upper triangular of comparisons to avoid duplicate
	// computation
	if(selected.isDiagonal()){
		//std::cerr << "diagnonal" << std::endl;
		this->n_comparisons_chunk = 0;
		for(U32 i = selected.fromRow; i < selected.toRow; ++i){
			for(U32 j = i; j < selected.toColumn; ++j){
				if(i == j) // upper triangular
					this->n_comparisons_chunk += (reader.getIndex().getContainer().at(i).n_variants * reader.getIndex().getContainer().at(i).n_variants - reader.getIndex().getContainer().at(i).n_variants) / 2;
				else // square
					this->n_comparisons_chunk += reader.getIndex().getContainer().at(i).n_variants * reader.getIndex().getContainer().at(j).n_variants;
			}
		}
		//std::cerr << "total variants chunk: " << this->n_comparisons_chunk << std::endl;
		//std::cerr << "parition size: " << this->n_comparisons_chunk / threads << std::endl;
		U64 n_variants_partition = this->n_comparisons_chunk / threads;
		if(threads == 1) n_variants_partition = std::numeric_limits<U64>::max();
		//std::cerr << "n_threads: " << threads << std::endl;

		U64 n_variants_counter = 0;
		std::pair<U32, U32> current;
		U32 currentThread = 0;

		for(U32 i = selected.fromRow; i < selected.toRow; ++i){
			current.first  = i;
			current.second = i;

			for(U32 j = i; j < selected.toColumn; ++j){
				if(i == j) n_variants_counter += (reader.getIndex().getContainer().at(i).n_variants * reader.getIndex().getContainer().at(i).n_variants - reader.getIndex().getContainer().at(i).n_variants) / 2;
				else n_variants_counter += reader.getIndex().getContainer().at(i).n_variants * reader.getIndex().getContainer().at(j).n_variants;

				if(n_variants_counter >= n_variants_partition){
					current.second = j + 1; // not inclusive range [A, B)
					if(current.second > selected.toColumn) current.second = selected.toColumn;
					//std::cerr << "break counter: " << current.first << "->" << current.second << " with " << n_variants_counter << "/" << n_variants_partition << " of " << this->n_comparisons_chunk << std::endl;
					this->thread_distribution[currentThread].push_back(LoadBalancerThread(i - selected.fromRow, current.first - selected.fromRow, current.second - selected.fromRow));

					current.first  = j + 1;
					current.second = j;
					n_variants_counter = 0;
					++currentThread;
				}
			}

			//std::cerr << current.second << std::endl;
			if(current.second != selected.toRow){
				current.second = selected.toRow;
				//std::cerr << "break new line: " << current.first << "->" << current.second << " with " << n_variants_counter <<"/" << n_variants_partition << std::endl;
				this->thread_distribution[currentThread].push_back(LoadBalancerThread(i - selected.fromRow, current.first - selected.fromRow, current.second - selected.fromRow));
				//parts.push_back(current);
				if(n_variants_counter >= n_variants_partition){
					//std::cerr << "break here" << std::endl;
					n_variants_counter = 0;
					++currentThread;
				}
			}
		}

		/*std::cerr << "currentThread: " << currentThread << std::endl;
		std::cerr << selected.toRow - selected.fromRow << std::endl;
		for(U32 i = 0; i < this->thread_distribution.size(); ++i){
			for(U32 j = 0; j < this->thread_distribution[i].size(); ++j){
				std::cerr << "idx: "<< i << "; row: " << this->thread_distribution[i][j].row << ":" << this->thread_distribution[i][j].fromColumn << "->" << this->thread_distribution[i][j].toColumn << "/" << reader.getIndex().getContainer().size() << std::endl;
			}
		}*/


	} else {
		//std::cerr << "square" << std::endl;
		const value_type& selectedRow = this->blocks[this->selected_chunk];
		const value_type& selectedCol = this->blocks[this->selected_chunk];

		this->n_comparisons_chunk = 0;
		for(U32 i = selectedRow.fromRow; i < selectedRow.toRow; ++i){
			for(U32 j = selectedCol.fromColumn; j < selectedCol.toColumn; ++j){
				this->n_comparisons_chunk += reader.getIndex().getContainer().at(i).n_variants * reader.getIndex().getContainer().at(j).n_variants;
			}
		}
		//std::cerr << "total variants chunk: " << this->n_comparisons_chunk << std::endl;
		//std::cerr << "parition size: " << this->n_comparisons_chunk / threads << std::endl;
		U64 n_variants_partition = this->n_comparisons_chunk / threads;
		if(threads == 1) n_variants_partition = std::numeric_limits<U64>::max();


		U64 n_variants_counter = 0;
		std::pair<U32, U32> current;
		U32 currentThread = 0;

		for(U32 i = selectedRow.fromRow; i < selectedRow.toRow; ++i){
			current.first  = selectedCol.fromColumn;
			current.second = selectedCol.fromColumn;

			for(U32 j = selectedCol.fromColumn; j < selectedCol.toColumn; ++j){
				n_variants_counter += reader.getIndex().getContainer().at(i).n_variants * reader.getIndex().getContainer().at(j).n_variants;

				if(n_variants_counter >= n_variants_partition){
					current.second = j + 1; // not inclusive range [A, B)
					if(current.second > selected.toColumn) current.second = selected.toColumn;
					//std::cerr << "break counter: " << current.first << "->" << current.second << " with " << n_variants_counter << "/" << n_variants_partition << std::endl;
					this->thread_distribution[currentThread].push_back(LoadBalancerThread(i - selectedRow.fromRow, (selectedRow.toColumn - selectedRow.fromColumn) + current.first - selectedCol.fromColumn, (selectedRow.toColumn - selectedRow.fromColumn) + current.second - selectedCol.fromColumn));

					current.first = j + 1;
					current.second = j;
					n_variants_counter = 0;
					++currentThread;
				}
			}

			//std::cerr << current.second << std::endl;
			if(current.second != selected.toColumn){
				current.second = selected.toColumn;
				//std::cerr << "break new line: " << current.first << "->" << current.second << std::endl;
				this->thread_distribution[currentThread].push_back(LoadBalancerThread(i - selectedRow.fromRow, (selectedRow.toColumn - selectedRow.fromColumn) + current.first - selectedCol.fromColumn, (selectedRow.toColumn - selectedRow.fromColumn) + current.second - selectedCol.fromColumn));
				//parts.push_back(current);
				if(n_variants_counter >= n_variants_partition){
					n_variants_counter = 0;
					++currentThread;
				}
			}
		}

		//std::cerr << "currentThread: " << currentThread << std::endl;
		//std::cerr << selected.toRow - selected.fromRow << std::endl;
		//for(U32 i = 0; i < this->thread_distribution.size(); ++i){
		//	for(U32 j = 0; j < this->thread_distribution[i].size(); ++j){
		//		std::cerr << "idx: "<< i << "; row: " << this->thread_distribution[i][j].row << ":" << this->thread_distribution[i][j].fromColumn << "->" << this->thread_distribution[i][j].toColumn << std::endl;
		//	}
		//}

	}

	//for(U32 i = 0; i < this->data_to_load.size(); ++i)
	//	std::cerr << i << '\t' << this->data_to_load[i].first << "->" << this->data_to_load[i].second << std::endl;

	return true;
}

bool LoadBalancerLD::setSelected(const S32 selected){
	if(selected < 0){
		std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Cannot set select a negative chunk..." << std::endl;
		return false;
	}

	this->selected_chunk = selected;
	return true;
}

bool LoadBalancerLD::setDesired(const S32 desired){
	if(desired < 0){
		std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Cannot cut workload into a negative number of blocks..." << std::endl;
		return false;
	}

	this->n_desired_chunks = desired;
	return true;
}

bool LoadBalancerLD::BuildInterval(const reader_type& reader, const U32 threads){
	if(reader.interval_tree_entries == nullptr){
		std::cerr << helpers::timestamp("ERR", "BALANCER") << "Cannot use slicing function without providing intervals..." << std::endl;
		return false;
	}

	// If slicing
	std::vector< std::pair<U32,U32> > target_blocks;

	for(U32 i = 0; i < reader.getHeader().getMagic().n_contigs; ++i){
		for(U32 j = 0; j < reader.interval_tree_entries[i].size(); ++j){
			std::cerr << reader.interval_tree_entries[i][j] << std::endl;
			std::pair<U32,U32> blocks = reader.getIndex().getContainer().findOverlap(reader.interval_tree_entries[i][j].contigID,
																	 reader.interval_tree_entries[i][j].start,
																	 reader.interval_tree_entries[i][j].stop);

			if(blocks.first == 0 && blocks.second == 0){
				std::cerr << helpers::timestamp("LOG", "BALANCER") << "No data found for interval (ID:" << reader.interval_tree_entries[i][j].contigID << ":" << reader.interval_tree_entries[i][j].start << "-" << reader.interval_tree_entries[i][j].stop << ")..." << std::endl;
				continue;
			}

			//std::cerr << "found blocks: " << blocks.first << "->" << blocks.second << std::endl;

			target_blocks.push_back(blocks);
			//for(U32 b = blocks.first; b < blocks.second; ++b){
			//	std::cerr << "Matches: " << b << std::endl;
			//}
		}
	}


	// Slicing does not overlap any data
	if(target_blocks.size() == 0){
		std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Could not find any data in the region(s) specified..." << std::endl;
		return false;
	}

	std::sort(target_blocks.begin(), target_blocks.end());
	std::vector< std::pair<U32,U32> > target_blocks_merged;
	target_blocks_merged.push_back(target_blocks[0]);
	std::cerr << target_blocks[0].first << "->" << target_blocks[0].second << std::endl;
	for(U32 i = 1; i < target_blocks.size(); ++i){
		std::cerr << target_blocks[i].first << "->" << target_blocks[i].second << std::endl;
		if(target_blocks[i].first <= target_blocks_merged.back().second)
			target_blocks_merged.back().second = target_blocks[i].second;
		else
			target_blocks_merged.push_back(target_blocks[i]);
	}

	std::cerr << "merged" << std::endl;
	for(U32 i = 0; i < target_blocks_merged.size(); ++i){
		std::cerr << target_blocks_merged[i].first << "->" << target_blocks_merged[i].second << std::endl;
	}

	// Count available blocks;
	U32 n_total_blocks = 0;
	for(U32 j = 0; j < target_blocks_merged.size(); ++j)
		n_total_blocks += target_blocks_merged[j].second - target_blocks_merged[j].first;

	std::cerr << "total blocks: " << n_total_blocks << std::endl;

	if(this->n_desired_chunks != 1){
		U32 cutSize = 1;
		for(U32 i = 1; i < n_total_blocks; ++i){
			if((i*i - i) / 2 + i == this->n_desired_chunks) // N choose 2 + N (upper triangular + diagonal)
				cutSize = i;
		}

		std::cerr << "n choose 2: n = " << cutSize << std::endl;
		if(cutSize == 1){
			std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Cannot cut into " << this->n_desired_chunks << " chunks. Chunks Have to be in the set choose(chunks,2) + chunks..." << std::endl;
			return(false);
		}

		const U32 blocks_per_partition = n_total_blocks / cutSize;
		std::cerr << "blocks/partition: " << blocks_per_partition << std::endl;
		U32 total = 0; // Sanity
		for(U32 i = 0; i < cutSize; ++i){ // column 0..N
			for(U32 j = i; j < cutSize; ++j){ // row c..N
				if(j + 1 == cutSize){
					if(i + 1 == cutSize){
						std::cerr << this->blocks.size() << ", last one: " << blocks_per_partition*i << "->" << n_total_blocks << ", " << blocks_per_partition*j << "->" << n_total_blocks << std::endl;
						this->blocks.push_back(value_type(blocks_per_partition*i, n_total_blocks, blocks_per_partition*j, n_total_blocks));
					} else {
						std::cerr << this->blocks.size() << ", last one: " << blocks_per_partition*i << "->" << blocks_per_partition*(i+1) << ", " << blocks_per_partition*j << "->" << n_total_blocks << std::endl;
						this->blocks.push_back(value_type(blocks_per_partition*i, blocks_per_partition*(i+1), blocks_per_partition*j, n_total_blocks));
					}
				} else {
					std::cerr << this->blocks.size() << ", normal: " << blocks_per_partition*i << "->" << blocks_per_partition*(i+1) << ", " << blocks_per_partition*j << "->" << blocks_per_partition*(j+1) << std::endl;
					this->blocks.push_back(value_type(blocks_per_partition*i, blocks_per_partition*(i+1), blocks_per_partition*j, blocks_per_partition*(j+1)));
				}
				++total;
			}
		}

		std::cerr << "First block: " << this->blocks[0] << std::endl;

		if(total != this->n_desired_chunks){
			std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Corrupted balancing..." << std::endl;
			return(false);
		}
	} else {
		this->blocks.push_back(value_type(0, n_total_blocks, 0, n_total_blocks));
	}



	// Data to load
	this->getSelectedLoad();
	for(U32 i = 0; i < this->data_to_load.size(); ++i){
		std::cerr << this->data_to_load[i].first << "->" << this->data_to_load[i].second << std::endl;
	}

	if(this->getSelectedLoadThreadsInterval(reader, threads, target_blocks_merged) == false){
		std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Corrupted balancing..." << std::endl;
		return(false);
	}

	std::cerr << "here after threads" << std::endl;

	if(this->getSelectedLoadInterval(target_blocks_merged) == false){
		std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Corrupted balancing..." << std::endl;
		return(false);
	}
	// Load absolute

	return(true);
}

bool LoadBalancerLD::Build(const reader_type& reader, const U32 threads){
	if(this->selected_chunk > this->n_desired_chunks){
		std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Incorrectly selected block (" << this->selected_chunk << '/' << this->n_desired_chunks << ")..." << std::endl;
		return false;
	}

	// Have intervals set
	if(reader.interval_tree_entries != nullptr)
		return(this->BuildInterval(reader, threads));

	// If selecting > 1 chunk
	if(this->n_desired_chunks != 1){
		U32 cutSize = 1;
		for(U32 i = 1; i < reader.getIndex().getContainer().size(); ++i){
			if((i*i - i) / 2 + i == this->n_desired_chunks) // N choose 2 + N (upper triangular + diagonal)
				cutSize = i;
		}

		if(cutSize == 1){
			std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Cannot cut into " << this->n_desired_chunks << " chunks. Chunks Have to be in the set choose(chunks,2) + chunks..." << std::endl;
			return(false);
		}

		// temp
		//std::cerr << helpers::timestamp("DEBUG") << "Cutsize is: " << cutSize << " and total blocks: " << reader.getIndex().getContainer().size() << std::endl;
		const U32 rowLength = reader.getIndex().getContainer().size() / cutSize;
		U32 total = 0; // Sanity
		for(U32 i = 0; i < cutSize; ++i){
			for(U32 j = i; j < cutSize; ++j){
				if(j + 1 == cutSize){
					//std::cerr << "last one: " << rowLength*i << "->" << reader.getIndex().getContainer().size() << ", " << rowLength*j << "->" << reader.getIndex().getContainer().size() << std::endl;
					this->blocks.push_back(value_type(rowLength*i, reader.getIndex().getContainer().size(), rowLength*j, reader.getIndex().getContainer().size()));
				} else {
					//std::cerr << "normal: " << rowLength*i << "->" << rowLength*(i+1) << ", " << rowLength*j << "->" << rowLength*(j+1) << std::endl;
					this->blocks.push_back(value_type(rowLength*i, rowLength*(i+1), rowLength*j, rowLength*(j+1)));
				}
				++total;
			}
		}

		if(total != this->n_desired_chunks){
			std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Corrupted balancing..." << std::endl;
			return(false);
		}

	} else {
		// All blocks
		this->blocks.push_back(value_type(0, reader.getIndex().getContainer().size(), 0, reader.getIndex().getContainer().size()));
	}

	// Data to load
	this->getSelectedLoad();

	// Get thread load
	if(!this->getSelectedLoadThreads(reader, threads))
		return false;

	return true;
}

bool LoadBalancerLD::BuildWindow(const reader_type& reader, const U32 threads, const U64 n_window_bases){
	if(this->selected_chunk > this->n_desired_chunks){
		std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Incorrectly selected block (" << this->selected_chunk << '/' << this->n_desired_chunks << ")..." << std::endl;
		return false;
	}

	if(n_window_bases < 1000){
		std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Probably not useful to set window size < 1000 bp" << std::endl;
		return false;
	}

	// temp
	U32 n_block_extension = 0;

	// If selecting > 1 chunk
	if(this->n_desired_chunks != 1){
		const U32 n_blocks_loaded = reader.getIndex().getContainer().size() / this->n_desired_chunks;
		const U32 from_block = n_blocks_loaded * this->selected_chunk;
		U32 to_block = n_blocks_loaded * (this->selected_chunk + 1);
		if(this->selected_chunk + 1 == this->n_desired_chunks) to_block = reader.getIndex().getContainer().size();

		// for i = 0:n_blocks_loaded
		// for j if not break
		//std::cerr << "Range is: " << from_block << "->" << to_block << "/" << reader.getIndex().getContainer().size() << std::endl;

		// Todo: extend command
		// Search for extension end
		if(this->selected_chunk + 1 != this->n_desired_chunks){
			const totempole::IndexEntry& index_entry_i = reader.getIndex().getContainer().at(to_block - 1);
			for(U32 j = to_block; j < reader.getIndex().getContainer().size(); ++j){
				const totempole::IndexEntry& index_entry_j = reader.getIndex().getContainer().at(j);
				if(index_entry_j.min_position - index_entry_i.max_position < n_window_bases){
					//std::cerr << "extend: " << j << " " << index_entry_j.min_position << "-" << index_entry_i.max_position << "=" << index_entry_j.min_position - index_entry_i.max_position << std::endl;
					++n_block_extension;
				} else
					break;
			}
		}

		this->data_to_load.push_back(std::pair<U32,U32>(from_block, to_block + n_block_extension));

	} else {
		// All blocks
		this->data_to_load.push_back(std::pair<U32,U32>(0, reader.getIndex().getContainer().size()));
	}

	//std::cerr << "Data to load: " << this->data_to_load[0].first << "->" << this->data_to_load[0].second << std::endl;
	//std::cerr << this->n_desired_chunks << std::endl;

	const U32 n_blocks = this->data_to_load[0].second - this->data_to_load[0].first + 1 - n_block_extension;
	U64 n_max_possible_comparisons = 0;
	for(U32 i = 0; i < n_blocks; ++i){
		const totempole::IndexEntry& index_entry_i = reader.getIndex().getContainer().at(i);
		n_max_possible_comparisons += (index_entry_i.n_variants*index_entry_i.n_variants + index_entry_i.n_variants) / 2 - index_entry_i.n_variants;
		//std::cerr << "keep internal" << std::endl;

		for(U32 j = i + 1; j < n_blocks + n_block_extension; ++j){
			const totempole::IndexEntry& index_entry_j = reader.getIndex().getContainer().at(j);

			// Check if they share the same contig
			if(index_entry_i.contigID != index_entry_j.contigID)
				break;

			if(index_entry_j.min_position - index_entry_i.max_position < n_window_bases){
				n_max_possible_comparisons += index_entry_j.n_variants * index_entry_i.n_variants;

				//std::cerr << "keep: " << i << "/" << j << " " << index_entry_j.min_position << "->" << index_entry_i.max_position << "=" << index_entry_j.min_position - index_entry_i.max_position << std::endl;
			} else {
				//std::cerr << "drop: " << i << "/" << j << " " << index_entry_j.min_position << "-" << index_entry_i.max_position << "=" << index_entry_j.min_position - index_entry_i.max_position << std::endl;
				break;
			}
		}
	}

	const U64 n_variants_thread = n_max_possible_comparisons / threads;
	//std::cerr << "Max total: " << n_max_possible_comparisons << "per thread " << n_variants_thread << std::endl;
	this->n_comparisons_chunk = n_max_possible_comparisons;

	this->thread_distribution.resize(threads);
	U64 n_total_current_thread = 0;
	U32 current_thread_id = 0;

	for(U32 i = 0; i < n_blocks; ++i){
		const totempole::IndexEntry& index_entry_i = reader.getIndex().getContainer().at(i);
		n_total_current_thread += (index_entry_i.n_variants*index_entry_i.n_variants + index_entry_i.n_variants) / 2 - index_entry_i.n_variants;

		U32 last_index = i;
		if(n_total_current_thread >= n_variants_thread){
			this->thread_distribution[current_thread_id].push_back(LoadBalancerThread(i, last_index, last_index + 1));
			n_total_current_thread = 0;
			++current_thread_id;
			last_index = i + 1;
		}

		bool breaking = false;

		for(U32 j = i + 1; j < n_blocks + n_block_extension; ++j){
			const totempole::IndexEntry& index_entry_j = reader.getIndex().getContainer().at(j);
			// Check if they share the same contig
			if(index_entry_i.contigID != index_entry_j.contigID)
				break;

			if(index_entry_j.min_position - index_entry_i.max_position < n_window_bases){
				n_total_current_thread += index_entry_j.n_variants * index_entry_i.n_variants;
			} else breaking = true;

			if(breaking || n_total_current_thread >= n_variants_thread){
				if(n_total_current_thread >= n_variants_thread){
					n_total_current_thread = 0;
					++current_thread_id;
				}
				//std::cerr << "Add: " << current_thread_id << ": " << i << "," << last_index << "->" << j + 1 << "/" << n_blocks << std::endl;
				this->thread_distribution[current_thread_id].push_back(LoadBalancerThread(i, last_index, j + !breaking));
				break;
			}
		}

		// Push back
		if(n_total_current_thread >= n_variants_thread && breaking == false){
			this->thread_distribution[current_thread_id].push_back(LoadBalancerThread(i, last_index, n_blocks));
			n_total_current_thread = 0;
			++current_thread_id;
		}
	}
	//std::cerr << current_thread_id << std::endl;
	//std::cerr << "total blocks: " << reader.getIndex().getContainer().size() << std::endl;

	//exit(1);

	// Get thread load
	//if(!this->getSelectedLoadThreads(reader, threads))
	//	return false;

	return true;
}

}
