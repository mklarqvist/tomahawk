#include "load_balancer_ld.h"

namespace Tomahawk{

LoadBalancerLD::LoadBalancerLD() : selected_chunk(0), desired_chunks(1){}
LoadBalancerLD::~LoadBalancerLD(){}

bool LoadBalancerLD::getSelectedLoad(){
	value_type selected = this->blocks[this->selected_chunk];

	// If there are both equal
	if(selected.fromRow == selected.fromColumn && selected.toRow == selected.toColumn){
		this->data_to_load.push_back(std::pair<U32, U32>(selected.fromRow, selected.toRow));
	} else {
		// No overlap
		this->data_to_load.push_back(std::pair<U32, U32>(selected.fromRow, selected.toRow));
		this->data_to_load.push_back(std::pair<U32, U32>(selected.fromColumn, selected.toColumn));
	}

	return true;
}

bool LoadBalancerLD::getSelectedLoadThreads(const reader_type& reader, const U32 threads){
	const value_type& selected = this->blocks[this->selected_chunk];
	this->thread_distribution.resize(threads);

	// limit
	if(selected.isDiagonal()){
		//std::cerr << "diagnonal" << std::endl;
		this->n_comparisons_chunk = 0;
		for(U32 i = selected.fromRow; i < selected.toRow; ++i){
			for(U32 j = i; j < selected.toColumn; ++j){
				if(i == j)
					this->n_comparisons_chunk += (reader.getIndex().getContainer().at(i).n_variants * reader.getIndex().getContainer().at(i).n_variants - reader.getIndex().getContainer().at(i).n_variants) / 2;
				else
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
		std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Cannot set select a negative chunk..." << std::endl;
		return false;
	}

	this->selected_chunk = selected;
	return true;
}

bool LoadBalancerLD::setDesired(const S32 desired){
	if(desired < 0){
		std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Cannot cut workload into a negative number of blocks..." << std::endl;
		return false;
	}

	this->desired_chunks = desired;
	return true;
}

bool LoadBalancerLD::Build(const reader_type& reader, const U32 threads){
	if(this->selected_chunk > this->desired_chunks){
		std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Incorrectly selected block (" << this->selected_chunk << '/' << this->desired_chunks << ")..." << std::endl;
		return false;
	}

	// If selecting > 1 chunk
	if(this->desired_chunks != 1){
		U32 cutSize = 1;
		for(U32 i = 1; i < reader.getIndex().getContainer().size(); ++i){
			if((i*i - i) / 2 + i == this->desired_chunks) // N choose 2 + N (upper triangular + diagonal)
				cutSize = i;
		}

		if(cutSize == 1){
			std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Cannot cut into " << this->desired_chunks << " chunks. Chunks Have to be in the set choose(chunks,2) + chunks..." << std::endl;
			return(false);
		}

		// temp
		//std::cerr << Helpers::timestamp("DEBUG") << "Cutsize is: " << cutSize << " and total blocks: " << reader.getIndex().getContainer().size() << std::endl;
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

		if(total != this->desired_chunks){
			std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Corrupted balancing..." << std::endl;
			return(false);
		}

	} else {
		// All blocks
		this->blocks.push_back(value_type(0, reader.getIndex().getContainer().size(), 0, reader.getIndex().getContainer().size()));
	}

	// Data to loadt
	this->getSelectedLoad();

	// Get thread load
	if(!this->getSelectedLoadThreads(reader, threads))
		return false;

	return true;
}

bool LoadBalancerLD::BuildWindow(const reader_type& reader, const U32 threads){
	if(this->selected_chunk > this->desired_chunks){
		std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Incorrectly selected block (" << this->selected_chunk << '/' << this->desired_chunks << ")..." << std::endl;
		return false;
	}

	// temp
	this->blocks.push_back(value_type(0, reader.getIndex().getContainer().size(), 0, reader.getIndex().getContainer().size()));
	//


	for(U32 i = 0; i < reader.getIndex().size(); ++i){
		const Totempole::IndexEntry& index_entry = reader.getIndex().getContainer().at(i);
		std::cerr << index_entry.contigID << "\t" << index_entry.min_position << "-" << index_entry.max_position << std::endl;
	}

	const U32 n_blocks = reader.getIndex().getContainer().size();
	U64 n_max_possible_comparisons = 0;
	for(U32 i = 0; i < n_blocks; ++i){
		const Totempole::IndexEntry& index_entry_i = reader.getIndex().getContainer().at(i);
		n_max_possible_comparisons += (index_entry_i.n_variants*index_entry_i.n_variants + index_entry_i.n_variants) / 2 - index_entry_i.n_variants;

		for(U32 j = i + 1; j < n_blocks; ++j){
			const Totempole::IndexEntry& index_entry_j = reader.getIndex().getContainer().at(j);
			if(index_entry_j.min_position - index_entry_i.max_position < 1e6){
				n_max_possible_comparisons += index_entry_j.n_variants * index_entry_i.n_variants;

				//std::cerr << "keep: " << i << "/" << j << " " << index_entry_j.min_position << "->" << index_entry_i.max_position << std::endl;
			} else {
				std::cerr << "drop: " << i << "/" << j << " " << index_entry_j.min_position << "-" << index_entry_i.max_position << "=" << index_entry_j.min_position - index_entry_i.max_position << std::endl;
				break;
			}
		}
	}

	const U64 n_variants_thread = n_max_possible_comparisons / threads;
	std::cerr << "Max total: " << n_max_possible_comparisons << "per thread " << n_variants_thread << std::endl;
	this->n_comparisons_chunk = n_max_possible_comparisons;

	this->thread_distribution.resize(threads);
	U64 n_total_current_thread = 0;
	U32 current_thread_id = 0;

	for(U32 i = 0; i < n_blocks; ++i){
		const Totempole::IndexEntry& index_entry_i = reader.getIndex().getContainer().at(i);
		n_total_current_thread += (index_entry_i.n_variants*index_entry_i.n_variants + index_entry_i.n_variants) / 2 - index_entry_i.n_variants;

		U32 last_index = i;
		if(n_total_current_thread >= n_variants_thread){
			this->thread_distribution[current_thread_id].push_back(LoadBalancerThread(i, last_index, last_index + 1));
			n_total_current_thread = 0;
			++current_thread_id;
			last_index = i + 1;
		}

		bool breaking = false;

		for(U32 j = i + 1; j < n_blocks; ++j){
			const Totempole::IndexEntry& index_entry_j = reader.getIndex().getContainer().at(j);
			if(index_entry_j.min_position - index_entry_i.max_position < 1e6){
				n_total_current_thread += index_entry_j.n_variants * index_entry_i.n_variants;
			} else {
				breaking = true;
			}

			if(breaking || n_total_current_thread >= n_variants_thread){
				if(n_total_current_thread >= n_variants_thread){
					n_total_current_thread = 0;
					++current_thread_id;
				}
				std::cerr << "Add: " << current_thread_id << ": " << i << "," << last_index << "->" << j + 1 << "/" << n_blocks << std::endl;
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
	std::cerr << current_thread_id << std::endl;
	std::cerr << "total blocks: " << reader.getIndex().getContainer().size() << std::endl;

	//exit(1);

	// Divide data into threads
	this->getSelectedLoad();

	// Get thread load
	//if(!this->getSelectedLoadThreads(reader, threads))
	//	return false;

	return true;
}

}
