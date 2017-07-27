#include <algorithm>
#include <bitset>
#include <queue>

#include "../../helpers.h"
#include "../../tomahawk/MagicConstants.h"
#include "TomahawkOutputReader.h"
#include "../../algorithm/sort/TomahawkOutputSort.h"

namespace Tomahawk {
namespace IO{
TomahawkOutputReader::TomahawkOutputReader() : position(0), size(0), interval_tree(nullptr), contigs(nullptr), contig_htable(nullptr), writer(nullptr), interval_tree_entries(nullptr)
{}

/*
std::cout << (int)(*this)[i].FLAGS << '\t' << (*this)[i].MAFMix << '\t' << (*this)[i].AcontigID << '\t' << (*this)[i].Aposition << '\t' << (*this)[i].BcontigID << '\t' << (*this)[i].Bposition
							<< '\t' << (*this)[i].p1 << '\t' << (*this)[i].p2 << '\t' << (*this)[i].q1 << '\t' << (*this)[i].q2 << '\t' << (*this)[i].D << '\t' << (*this)[i].Dprime
							<< '\t' << (*this)[i].R2 << '\t' << (*this)[i].P << '\t' << (*this)[i].chiSqFisher << '\t' << (*this)[i].chiSqModel << '\n';
*/


bool TomahawkOutputReader::view(const std::string& input){
	if(this->interval_tree != nullptr)
		return(this->__viewRegion());
	else
		return(this->__viewFilter());
}

bool TomahawkOutputReader::__viewRegion(void){
	if(this->interval_tree != nullptr){
		const Tomahawk::IO::TomahawkOutputEntry*  entry;
		std::vector<interval_type> rets;
		while(this->nextVariant(entry)){
			// If iTree for contigA exists
			if(this->interval_tree[entry->AcontigID] != nullptr){
				rets = this->interval_tree[entry->AcontigID]->findOverlapping(entry->Aposition, entry->Aposition);
				if(rets.size() > 0){
					for(U32 i = 0; i < rets.size(); ++i){
						if(rets[i].value != nullptr){
							if((entry->BcontigID == rets[i].value->contigID) &&
							   (entry->Bposition >= rets[i].value->start && entry->Bposition <= rets[i].value->stop)){
								//std::cerr << "hit linked A" << std::endl;
								if(this->filter.filter(*entry))
									entry->write(std::cout, this->contigs);

								break;
							}
						} else {
							if(this->filter.filter(*entry))
								entry->write(std::cout, this->contigs);

							break;
						}
					}
					continue;
				}
			}

			// If iTree for contigB exists
			if(this->interval_tree[entry->BcontigID] != nullptr){
				rets = this->interval_tree[entry->BcontigID]->findOverlapping(entry->Bposition, entry->Bposition);
				if(rets.size() > 0){
					for(U32 i = 0; i < rets.size(); ++i){
						if(rets[i].value != nullptr){
							if((entry->AcontigID == rets[i].value->contigID) &&
							   (entry->Aposition >= rets[i].value->start && entry->Aposition <= rets[i].value->stop)){
								//std::cerr << "hit linked B" << std::endl;
								if(this->filter.filter(*entry))
									entry->write(std::cout, this->contigs);

								break;
							}
						} else {
							if(this->filter.filter(*entry))
								entry->write(std::cout, this->contigs);

							break;
						}
					}
					continue;
				}
			}
		}
	}

	return true;
}

bool TomahawkOutputReader::__viewOnly(void){
	std::cerr << "illegal" << std::endl;
	exit(1);

	const entry_type* entry;
	while(this->reader.nextEntry(entry)){
		std::cout << *entry << '\n';
		//std::cout.write(reinterpret_cast<const char*>(entry), sizeof(entry_type));
	}

	return true;
}

bool TomahawkOutputReader::__viewFilter(void){
	const Tomahawk::IO::TomahawkOutputEntry*  entry;
	while(this->nextVariant(entry)){
		if(this->filter.filter(*entry))
//			std::cout << this->contigs[entry->AcontigID].name << '\t' << this->contigs[entry->BcontigID].name << '\t' << *entry << '\n';
			entry->write(std::cout, this->contigs);
	}

	return true;
}

bool TomahawkOutputReader::AddRegions(std::vector<std::string>& positions){
	if(positions.size() == 0)
		return false;

	if(this->interval_tree == nullptr){
		this->interval_tree = new tree_type*[this->header.n_contig];
		for(U32 i = 0; i < this->header.n_contig; ++i)
			this->interval_tree[i] = nullptr;
	}
	if(this->interval_tree_entries == nullptr)
		this->interval_tree_entries = new std::vector<interval_type>[this->header.n_contig];

	for(U32 i = 0; i < positions.size(); ++i){
		//std::cerr << i << ": " << positions[i] << std::endl;
		// Pattern cA:pAf-pAt;cB:pBf-pBt
		if(positions[i].find(',') != std::string::npos){
			//std::cerr << "linked intervals" << std::endl;
			std::vector<std::string> ret = Helpers::split(positions[i], ',');
			if(ret.size() == 1){
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
				return false;

			} else if(ret.size() == 2){
				// parse left
				interval_type intervalLeft;
				if(!__ParseRegion(ret[0], intervalLeft))
					return false;

				// parse right
				interval_type intervalRight;
				if(!__ParseRegion(ret[1], intervalRight))
					return false;

				// Todo: WARNING
				// This results in illegal pointers if the vector resizes
				// and pointers change
				this->interval_tree_entries[intervalLeft.contigID].push_back(interval_type(intervalLeft));
				this->interval_tree_entries[intervalRight.contigID].push_back(interval_type(intervalRight));
				if(intervalLeft.contigID != intervalRight.contigID){
					this->interval_tree_entries[intervalLeft.contigID].back().value = &this->interval_tree_entries[intervalRight.contigID].back();
					this->interval_tree_entries[intervalRight.contigID].back().value = &this->interval_tree_entries[intervalLeft.contigID].back();

					// Link the intervals together
					//std::cerr << this->interval_tree_entries[intervalLeft.contigID].back() << '\t' << this->interval_tree_entries[intervalRight.contigID].back() << std::endl;
				} else {
					this->interval_tree_entries[intervalLeft.contigID].back().value = &this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2];
					this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2].value = &this->interval_tree_entries[intervalLeft.contigID].back();

					//std::cerr << this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2] << '\t' << this->interval_tree_entries[intervalRight.contigID].back() << std::endl;
				}

			} else {
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
				return false;
			}
		} else {
			interval_type interval;
			if(!__ParseRegion(positions[i], interval))
				return false;

			this->interval_tree_entries[interval.contigID].push_back(interval_type(interval));
		}
	}

	for(U32 i = 0; i < this->header.n_contig; ++i){
		if(this->interval_tree_entries[i].size() != 0){
			//std::cerr << "constructing itree for id: " << i << std::endl;
			this->interval_tree[i] = new tree_type(this->interval_tree_entries[i]);
		} else
			this->interval_tree[i] = nullptr;
	}

	return true;
}

bool TomahawkOutputReader::__ParseRegion(const std::string& region, interval_type& interval){
	std::vector<std::string> ret = Helpers::split(region, ':');
	if(ret.size() == 1){
		if(ret[0].find('-') != std::string::npos){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}

		// is contigID only
		//std::cerr << "contigONly" << std::endl;
		U32* contigID;
		if(!this->contig_htable->GetItem(&region[0], &region, contigID, region.size())){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << region << " is not defined in the header!" << std::endl;
			return false;
		}
		interval(0, this->contigs[*contigID].bases, *contigID);

	} else if(ret.size() == 2){
		// is contigID:pos-pos
		U32* contigID;
		if(!this->contig_htable->GetItem(&ret[0][0], &ret[0], contigID, ret[0].size())){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << ret[0] << " is not defined in the header!" << std::endl;
			return false;
		}

		std::vector<std::string> retPos = Helpers::split(ret[1], '-');
		if(retPos.size() == 1){
			// only one pos
			const double pos = std::stod(retPos[0]);
			//std::cerr << "single position: " << pos << std::endl;
			interval(pos, pos, *contigID);

		} else if(retPos.size() == 2){
			// is two positions
			double posA = std::stod(retPos[0]);
			double posB = std::stod(retPos[1]);
			if(posB < posA){
				//std::cerr << "end position > start position: swapping" << std::endl;
				std::swap(posA, posB);
			}
			//std::cerr << "full region: " << this->contigs[*contigID].name << ":" << posA << '-' << posB << std::endl;
			interval(posA, posB, *contigID);

		} else {
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}
	} else {
		std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::Open(const std::string input){
	this->stream.open(input, std::ios::binary | std::ios::in | std::ios::ate);
	if(!this->stream.good()){
		std::cerr << "failed to open file " << input << std::endl;
		return false;
	}

	this->filesize = this->stream.tellg();
	this->stream.seekg(0);

	if(!this->stream.good()){
		std::cerr << "bad stream" << std::endl;
		return false;
	}

	this->stream >> this->header;
	if(!this->header.validate(Tomahawk::Constants::WRITE_HEADER_LD_MAGIC)){
		std::cerr << "failed to validate header" << std::endl;
		return false;
	}

	if(!this->ParseHeader()){
		std::cerr << "failed to parse header" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::ParseHeader(void){
	if(this->header.n_contig == 0)
		return false;

	if(this->header.n_contig < 1024)
		this->contig_htable = new hash_table(1024);
	else
		this->contig_htable = new hash_table(this->header.n_contig * 2);

	this->contigs = new contig_type[this->header.n_contig];
	U32* ret;
	for(U32 i = 0; i < this->header.n_contig; ++i){
		this->stream >> this->contigs[i];
		// std::cerr << this->contigs[i] << std::endl;
		if(!this->contig_htable->GetItem(&this->contigs[i].name[0], &this->contigs[i].name, ret, this->contigs[i].name.size())){
			// Add to hash table
			this->contig_htable->SetItem(&this->contigs[i].name[0], &this->contigs[i].name, i, this->contigs[i].name.size());
		} else {
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Duplicated contig name: " << this->contigs[i].name << "!" << std::endl;
			exit(1); // unrecoverable error
		}
	}

	return true;
}

bool TomahawkOutputReader::nextBlock(void){
	// Stream died
	if(!this->stream.good()){
		std::cerr << "stream died" << std::endl;
		return false;
	}

	// EOF
	if(this->stream.tellg() == this->filesize){
		//std::cerr << "eof" << std::endl;
		return false;
	}

	buffer.resize(sizeof(tgzf_type));
	this->stream.read(&buffer.data[0],  Constants::TGZF_BLOCK_HEADER_LENGTH);
	const tgzf_type* h = reinterpret_cast<const tgzf_type*>(&buffer.data[0]);
	buffer.pointer = Constants::TGZF_BLOCK_HEADER_LENGTH;
	if(!h->Validate()){
		std::cerr << "failed to validate" << std::endl;
		return false;
	}

	buffer.resize(h->BSIZE);

	// Recast because if buffer is resized then the pointer address is incorrect
	// resulting in segfault
	h = reinterpret_cast<const tgzf_type*>(&buffer.data[0]);

	this->stream.read(&buffer.data[Constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - Constants::TGZF_BLOCK_HEADER_LENGTH);
	if(!this->stream.good()){
		std::cerr << "truncated file" << std::endl;
		return false;
	}

	buffer.pointer = h->BSIZE;
	const U32* outsize = reinterpret_cast<const U32*>(&buffer[buffer.pointer -  sizeof(U32)]);
	//const U32* crc = reinterpret_cast<const U32*>(&buffer[buffer.pointer -  sizeof(U32) - sizeof(U32)]);
	//std::cerr << *outsize << '\t' << *crc << std::endl;
	output_buffer.resize(*outsize);
	this->output_buffer.reset();

	if(!this->gzip_controller.Inflate(buffer, output_buffer)){
		std::cerr << "failed inflate" << std::endl;
		return false;
	}

	if(this->output_buffer.size() == 0){
		std::cerr << "empty data" << std::endl;
		return false;
	}

	// Reset buffer
	this->buffer.reset();

	// Reset iterator position and size
	this->position = 0;
	this->size = this->output_buffer.size() / sizeof(TomahawkOutputEntry);

	// Validity check
	if(this->output_buffer.size() % sizeof(TomahawkOutputEntry) != 0){
		std::cerr << "data is corrupted" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::nextVariant(const TomahawkOutputEntry*& entry){
	if(this->position == this->size){
		if(!this->nextBlock())
			return false;
	}

	entry = (*this)[this->position];
	++this->position;

	return true;
}

bool TomahawkOutputReader::summary(const std::string& input){
	if(!this->reader.setup(input))
		return false;

	return true;
}

bool TomahawkOutputReader::index(const std::string& input){
	if(!this->reader.setup(input))
		return false;

	const entry_type* entry;
	if(!this->reader.nextEntry(entry))
		return false;

	//const entry_type* previous = entry;
	U32 currentAID = entry->AcontigID;
	U32 currentAPos = entry->Aposition;
	U32 currentBID = entry->BcontigID;
	U64 AIDSteps = 0;
	U64 APosSteps = 0;
	U64 BIDSteps = 0;

	double AposStepsR = 0;

	U64 outputEntries = 0;

	while(this->reader.nextEntry(entry)){
		/*
		if(!(*previous < *entry)){
			std::cerr << "file is not sorted" << std::endl;
			std::cerr << previous->AcontigID << '\t' << previous->Aposition << '\t' << previous->BcontigID << '\t' << previous->Bposition << std::endl;
			std::cerr << entry->AcontigID << '\t' << entry->Aposition << '\t' << entry->BcontigID << '\t' << entry->Bposition << std::endl;
			return false;
		}
		std::swap(entry, previous);
		*/
		++AIDSteps;
		++APosSteps;
		++BIDSteps;
		AposStepsR += entry->R2;

		if(entry->BcontigID != currentBID || entry->Aposition != currentAPos || entry->AcontigID != currentAID){
			//std::cerr << "switch: " << currentAID << ',' << currentAPos << ',' << currentBID << '\t' << entry->AcontigID << ',' << entry->Aposition << ',' << entry->BcontigID << '\t' << BIDSteps << std::endl;
			currentBID = entry->BcontigID;
			BIDSteps = 0;
			++outputEntries;
		}

		if(entry->Aposition != currentAPos || entry->AcontigID != currentAID){
			std::cout << currentAID << '\t' << currentAPos << '\t' << APosSteps << '\t' << AposStepsR/APosSteps << '\n';
			currentAPos = entry->Aposition;
			APosSteps = 0;
			AposStepsR = 0;
			++outputEntries;
		}

		if(entry->AcontigID != currentAID){
			//std::cerr << "switch: " << currentAID << "->" << entry->AcontigID << '\t' << AIDSteps << std::endl;
			currentAID = entry->AcontigID;
			AIDSteps = 0;
			++outputEntries;
		}

	}
	std::cerr << "Index would have: " << outputEntries << " entries..." << std::endl;

	return true;
}


}
} /* namespace Tomahawk */
