#include "TomahawkCalculations.h"

namespace Tomahawk {

TomahawkCalculations::TomahawkCalculations() : group_htable(nullptr)
{
}

TomahawkCalculations::~TomahawkCalculations()
{
	delete this->group_htable;
}

bool TomahawkCalculations::loadGroups(const std::string& file){
	if(this->Occ.size() != 0){
		std::cerr << "groups already loaded" << std::endl;
		return false;
	}

	if(file.size() == 0){
		std::cerr << "No file set" << std::endl;
		return false;
	}

	// Open stream
	std::ifstream stream(file);
	if(!stream.good()){
		std::cerr << "bad file" << std::endl;
		return false;
	}

	this->group_htable = new hash_table(2048);


	std::string line;
	U32 n_lines = 1;

	// Grab first line
	if(!getline(stream, line)){
		std::cerr << "failed to get first line" << std::endl;
		return false;
	}

	// count tabs in first line
	// Assert correct format
	// Count tabs until out-of-range
	size_t pos = -1;
	U32 tabs = 0;
	while(true){
		pos = line.find('\t', pos + 1);
		if(pos == std::string::npos)
			break;

		++tabs;
	}
	++tabs; // eof tab

	std::vector< std::vector< S32 > > groupings(this->samples);

	S32* sampleID = nullptr;
	S32* groupID_lookup = nullptr;
	S32 groupID = 0;

	while(getline(stream, line)){
		// Empty lines
		if(line.size() == 0)
			break;

		// Assert correct format
		// Count tabs until out-of-range
		size_t prev_pos = 0;
		U32 inner_tabs = 0;
		size_t pos = line.find('\t', prev_pos + 1);


		const std::string sampleName = std::string(&line[prev_pos], pos - prev_pos);

		if(!this->totempole.sampleHashTable->GetItem(&sampleName[0], &sampleName, sampleID, sampleName.length())){
			std::cerr << "sample does not exist" << std::endl;
			return false;
		}
		//std::cerr << sampleName << '\t' << *sampleID << '\t';

		prev_pos = pos + 1;
		++inner_tabs;

		for(U32 i = 1; i < tabs; ++i){
			size_t pos = line.find('\t', prev_pos + 1);
			if(pos == std::string::npos){
				if(i + 1 != tabs){
					std::cerr << "mangled data" << std::endl;
					return false;
				}
				pos = line.size();
			}
			//std::cerr << pos << '\t' << prev_pos << std::endl;
			const std::string group = std::string(&line[prev_pos], pos - prev_pos);
			if(!this->group_htable->GetItem(&group[0], &group, groupID_lookup, group.length())){
				this->group_htable->SetItem(&group[0], &group, groupID, group.length());
				this->groups.push_back(group);

				//std::cerr << "added: " << group << " with id " << groupID << std::endl;
				groupings[*sampleID].push_back(groupID);
				++groupID;
			} else {
				groupings[*sampleID].push_back(*groupID_lookup);
			}

			//std::cerr << std::string(&line[prev_pos], pos - prev_pos) << '\t';
			prev_pos = pos + 1;

			++inner_tabs;
		}

		if(tabs != inner_tabs){
			std::cerr << Helpers::timestamp("ERROR") << "Illegal format! Expected " << tabs << " columns! Line: " << n_lines << "..." << std::endl;
			return false;
		}

		//std::cerr << std::endl;

		// do stuff
		++n_lines;
	}

	if(groupID == 0){
		std::cerr << "no data loaded" << std::endl;
		return false;
	}

	this->Occ = occ_matrix(this->samples + 1, std::vector< U64 >(groupID, 0));
	occ_vector* prev = &Occ[0];

	for(U32 i = 0; i < this->samples; ++i){
		// Propagate previous vector counts
		for(U32 j = 0; j < groupID; ++j)
			this->Occ[i + 1][j] = prev->at(j);

		// Update
		for(U32 j = 0; j < groupings[i].size(); ++j)
			this->Occ[i + 1][groupings[i][j]] = prev->at(groupings[i][j]) + 1;

		prev = &this->Occ[i + 1];
	}

	// Temp
	// Dump
	/*
	for(U32 i = 0; i < this->Occ.size(); ++i){
		for(U32 j = 0; j < this->Occ[0].size(); ++j){
			std::cout << this->Occ[i][j] << '\t';
		}
		std::cout << std::endl;
	}
	*/

	return true;
}


bool TomahawkCalculations::calculateTajimaD(const U32 bin_size){
	switch(this->bit_width){
	case 1: return(this->__calculateTajimaD<BYTE>(bin_size));
	case 2: return(this->__calculateTajimaD<U16>(bin_size));
	case 4: return(this->__calculateTajimaD<U32>(bin_size));
	case 8: return(this->__calculateTajimaD<U64>(bin_size));
	default: exit(1); break;
	}

	return false;
}

bool TomahawkCalculations::calculateFST(void){
	switch(this->bit_width){
	case 1: return(this->__calculateFST<BYTE>());
	case 2: return(this->__calculateFST<U16>());
	case 4: return(this->__calculateFST<U32>());
	case 8: return(this->__calculateFST<U64>());
	default: exit(1); break;
	}

	return false;
}

bool TomahawkCalculations::calculateSFS(void){
	switch(this->bit_width){
	case 1: return(this->__calculateSFSGrouped<BYTE>());
	case 2: return(this->__calculateSFSGrouped<U16>());
	case 4: return(this->__calculateSFSGrouped<U32>());
	case 8: return(this->__calculateSFSGrouped<U64>());
	default: exit(1); break;
	}

	return false;
}

} /* namespace Tomahawk */
