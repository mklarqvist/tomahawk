#ifndef ALGORITHM_BALANCER_H_
#define ALGORITHM_BALANCER_H_

#include "LoadBalancerBlock.h"

namespace Tomahawk{

class Balancer{
	typedef LoadBalancerBlock block_type;

public:
	Balancer() : selected_chunk(0), desired_chunks(1){}
	~Balancer(){}

	bool getSelectedLoad(){
		//std::cerr << Helpers::timestamp("DEBUG", "BALANCER") << "What data to load?" << std::endl;
		block_type selected = this->blocks[this->selected_chunk];
		//std::cerr << Helpers::timestamp("DEBUG", "BALANCER") << this->selected_chunk << '/' << this->blocks.size() << std::endl;

		// attempt to merge
		// If there are both equal
		if(selected.fromRow == selected.fromColumn && selected.toRow == selected.toColumn){
			this->data_to_load.push_back(std::pair<U32, U32>(selected.fromRow, selected.toRow));
			//std::cerr << "same: " << selected << std::endl;
		} else {
			// No voerlap
			//std::cerr << "Not same: " << selected << std::endl;
			this->data_to_load.push_back(std::pair<U32, U32>(selected.fromRow, selected.toRow));
			this->data_to_load.push_back(std::pair<U32, U32>(selected.fromColumn, selected.toColumn));
		}

		return true;
	}

	bool getSelectedLoadThreads(const U32 threads){
		const block_type& selected = this->blocks[selected_chunk];
		//std::cerr << Helpers::timestamp("DEBUG", "BALANCER") << "Thread balancing..." << std::endl;

		this->thread_distribution.resize(threads);

		if(threads == 1){
			this->thread_distribution[0].push_back(block_type(0, selected.getRows(), 0, selected.getColumns(), selected.fromRow, selected.toRow, selected.fromColumn, selected.toColumn, selected.isDiagonal()));
			return true;
		}

		//
		if(selected.isDiagonal()){
			if(!SILENT){
				std::cerr << Helpers::timestamp("LOG", "BALANCER") << "Case is diagonal (chunk " << this->selected_chunk << '/' << this->desired_chunks << ")..." << std::endl;
				std::cerr << Helpers::timestamp("LOG", "BALANCER") << "Total comparisons: " << selected.getSize() << " and per thread: " << selected.getSize()/threads << std::endl;
			}

			U32 loadThread = selected.getSize()/threads;
			U32 it = 0;
			U32 from = 0;
			U32 fromCol = 0;
			U32 threadID = 0;

			//
			for(U32 i = 0; i < selected.getRows(); ++i){
				for(U32 j = i; j < selected.getColumns(); ++j){
					++it;

					// If number of comparions over threshold
					if(it >= loadThread){
						// if broken over a line
						// i.e. not broken on the same line number
						if(from == i){
							//std::cerr << "B\t" << threadID << ": " << from << '-' << i+1 << '\t' << fromCol << '-' << j << '\t' << selected.fromRow+from << '-' << selected.fromRow+(i+1) << '\t' << selected.fromColumn+fromCol << '-' << selected.toColumn+j << std::endl;
							this->thread_distribution[threadID].push_back(block_type(from, i+1, fromCol, j, selected.fromRow+from, selected.fromRow+i+1, selected.fromColumn+fromCol, selected.fromColumn+j));
						}
						// If broken over multiple lines
						else {
							if(threadID + 1 == threads){
								i = selected.getRows() - 1;
								j = selected.getColumns();
							}

							// If next line: no middle full lines
							if(from + 1 == i){
								//std::cerr << "N\t" << threadID << ": " << from << '-' << from+1 << '\t' << fromCol << '-' << selected.getColumns() << '\t' << "FALSE" << std::endl;
								//std::cerr << "N\t" << threadID << ": " << i << '-' << i+1 << '\t' << i << '-' << j << '\t' << "FALSE" << std::endl;
								this->thread_distribution[threadID].push_back(block_type(from, from+1, fromCol, selected.getColumns(), selected.fromRow+from, selected.fromRow+from+1, selected.fromColumn+fromCol, selected.toColumn));
								this->thread_distribution[threadID].push_back(block_type(i, i+1, i, j, selected.fromRow+i, selected.fromRow+i+1, selected.fromColumn+i, selected.fromColumn+j));
								fromCol = j;
								from = i;
							} else {
								//std::cerr << "E\t" << threadID << ": " << from << '-' << from + 1 << '\t' << fromCol << '-' << selected.getColumns() << '\t' << selected.fromRow+from << '-' << selected.fromRow+(from+1) << '\t' << selected.fromColumn+fromCol << '-' << selected.toColumn << std::endl;
								//std::cerr << "E\t" << threadID << ": " << from + 1 << '-' << i << '\t' << from + 1 << '-' << selected.getColumns() << '\t' << selected.fromRow+from+1 << '-' << selected.fromRow+(i) << '\t' << selected.fromColumn+from+1 << '-' << selected.toColumn << std::endl;
								//std::cerr << "E\t" << threadID << ": " << i << '-' << i + 1 << '\t' << i << '-' << j << '\t' << selected.fromRow+i << '-' << selected.fromRow+(i+1) << '\t' << selected.fromColumn+i << '-' << selected.fromColumn+j << std::endl;
								this->thread_distribution[threadID].push_back(block_type(from, from + 1, fromCol, selected.getColumns(), selected.fromRow+from, selected.fromRow+from+1, selected.fromColumn+fromCol, selected.toColumn));
								this->thread_distribution[threadID].push_back(block_type(from + 1, i, from + 1, selected.getColumns(), selected.fromRow+from+1, selected.fromRow+i, selected.fromColumn+from+1, selected.toColumn, true));
								this->thread_distribution[threadID].push_back(block_type(i, i + 1, i, j, selected.fromRow+i, selected.fromRow+i+1, selected.fromColumn+i, selected.fromColumn+j));
							}
						}
						it = 0;
						from = i;
						fromCol = j;
						++threadID;
					}

				}
			}
		}
		// Is not a diagonal square
		else {
			if(!SILENT){
				std::cerr << Helpers::timestamp("LOG", "BALANCER") << "Case is square (chunk " << this->selected_chunk << '/' << this->desired_chunks << ")..." << std::endl;
				std::cerr << Helpers::timestamp("LOG", "BALANCER") << "Total comparisons: " << selected.getSize() << " and per thread: " << selected.getSize()/threads << std::endl;
			}

			U32 loadThread = selected.getSize()/threads;
			U32 it = 0;
			U32 from = 0;
			U32 fromCol = selected.getRows();
			U32 threadID = 0;

			//
			for(U32 i = 0; i < selected.getRows(); ++i){
				for(U32 j = selected.getRows(); j < 2*selected.getRows(); ++j){
					++it;

					// If number of comparions over threshold
					if(it >= loadThread){
						// if broken over a line
						// i.e. not broken on the same line number
						if(from == i){
							//std::cerr << threadID << ": " << from << '-' << i+1 << '\t' << fromCol << '-' << j << '\t' << "FALSE" << std::endl;
							this->thread_distribution[threadID].push_back(block_type(from, i+1, fromCol, j, selected.fromRow+from, selected.fromRow+i+1, selected.fromColumn+fromCol, selected.fromColumn+j));
						}
						// If broken over multiple lines
						else {
							if(threadID + 1 == threads){
								i = selected.getRows() - 1;
								j = 2*selected.getRows();
							}

							// If next line: no middle full lines
							if(from + 1 == i){
								//std::cerr << threadID << ": " << from << '-' << from+1 << '\t' << fromCol << '-' << 2*selected.getRows() << '\t' << "FALSE" << std::endl;
								//std::cerr << threadID << ": " << i << '-' << i+1 << '\t' << selected.getRows() << '-' << j << '\t' << "FALSE" << std::endl;
								this->thread_distribution[threadID].push_back(block_type(from, from+1, fromCol, 2*selected.getRows(), selected.fromRow+from, selected.fromRow+from+1, selected.fromColumn+fromCol, selected.fromColumn+2*selected.getRows()));
								this->thread_distribution[threadID].push_back(block_type(i, i+1, 2*selected.getRows(), j, selected.fromRow+i, selected.fromRow+i+1, selected.fromColumn+2*selected.getRows(), selected.fromColumn+j));
								fromCol = j;
								from = i;
							} else {
								//std::cerr << threadID << ": " << from << '-' << from + 1 << '\t' << fromCol << '-' << 2*selected.getRows() << '\t' << "FALSE" << std::endl;
								//std::cerr << threadID << ": " << from + 1 << '-' << i << '\t' << selected.getRows() << '-' << 2*selected.getRows() << '\t' << "FALSE" << std::endl;
								//std::cerr << threadID << ": " << i << '-' << i + 1 << '\t' << selected.getRows() << '-' << j << '\t' << "FALSE" << std::endl;
								this->thread_distribution[threadID].push_back(block_type(from, from + 1, fromCol, 2*selected.getRows(), selected.fromRow+from, selected.fromRow+from+1, selected.fromColumn+fromCol, selected.fromColumn+2*selected.getRows()));
								this->thread_distribution[threadID].push_back(block_type(from + 1, i, selected.getRows(), 2*selected.getRows(), selected.fromRow+from+1, selected.fromRow+i, selected.fromColumn+selected.getRows(), selected.fromColumn+2*selected.getRows()));
								this->thread_distribution[threadID].push_back(block_type(i, i + 1, selected.getRows(), j, selected.fromRow+i, selected.fromRow+i+1, selected.fromColumn+selected.getRows(), selected.fromColumn+j));
							}
						}
						it = 0;
						from = i;
						fromCol = j;
						++threadID;
					}

				}
			}
		}

		// assertion


		//std::cerr << "DEBUG" << std::endl;
		//for(U32 i = 0; i < this->thread_distribution.size(); ++i)
		//	std::cerr << i << '\t' << this->thread_distribution[i].size() << std::endl;
		//std::cerr << "Has: " << this->thread_distribution.size() << " thread blocks" << std::endl;

		return true;
	}

	bool setSelected(const S32 selected){
		if(selected < 0){
			std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Cannot set select a negative chunk..." << std::endl;
			return false;
		}

		this->selected_chunk = selected;
		return true;
	}

	bool setDesired(const S32 desired){
		if(desired < 0){
			std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Cannot cut workload into a negative number of blocks..." << std::endl;
			return false;
		}

		this->desired_chunks = desired;
		return true;
	}

	bool Build(const U32 total_blocks, const U32 threads){
		if(this->selected_chunk > this->desired_chunks){
			std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Incorrectly selected block (" << this->selected_chunk << '/' << this->desired_chunks << ")..." << std::endl;
			return false;
		}

		// If selecting > 1 chunk
		if(this->desired_chunks != 1){
			U32 cutSize = 1;
			//std::vector<U32> backup_cuts;
			for(U32 i = 1; i < total_blocks; ++i){

				if((i*i - i) / 2 == this->desired_chunks)
					cutSize = i;

			}

			if(cutSize == 1){
				std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Cannot cut into " << this->desired_chunks << " chunks" << std::endl;
				return(false);
			}

			U32 total = 0; // Sanity
			//std::cerr << "cut-size is: " << cutSize << std::endl;
			const U32 rowLength = total_blocks / cutSize;
			for(U32 i = 0; i < cutSize-1; ++i){
				//std::cerr << i << '/' << cutSize-1 << std::endl;
				U32 j = i;
				U32 fromX = i*rowLength;
				U32 toX = (i+1)*rowLength;
				if(i + 1 == cutSize - 1)
					toX = total_blocks;

				for(; j < cutSize-1; ++j){
					U32 fromY = j*rowLength;
					U32 toY = (j+1)*rowLength;


					if(j + 1 == cutSize - 1)
						toY = total_blocks;

					//std::cerr << "(" << i << ',' << j << ")\t" << fromX << '-' << toX << '\t' << fromY << '-' << toY << std::endl;
					this->blocks.push_back(block_type(fromX, toX, fromY, toY));
					++total;
				}
			}

			//std::cerr << "Total: " << total << '/' << this->desired_chunks << std::endl;
			if(total != this->desired_chunks){
				std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Corrupted balancing..." << std::endl;
				return(false);
			}

		} else {
			// All blocks
			this->blocks.push_back(block_type(0, total_blocks, 0, total_blocks));
		}

		// What data do we load?
		this->getSelectedLoad();

		// Divide data into threads
		if(!this->getSelectedLoadThreads(threads))
			return false;

		return true;
	}

	inline std::vector< std::pair<U32, U32> >& getLoad(void){ return(this->data_to_load); }

public:
	U32 selected_chunk;
	U32 desired_chunks;

	std::vector<block_type> blocks;
	std::vector< std::pair<U32, U32> > data_to_load;
	std::vector< std::vector<block_type> > thread_distribution;
};

}
#endif /* ALGORITHM_BALANCER_H_ */
