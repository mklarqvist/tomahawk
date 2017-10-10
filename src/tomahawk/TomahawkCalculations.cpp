#include "TomahawkCalculations.h"

namespace Tomahawk {

TomahawkCalculations::TomahawkCalculations()
{
}

TomahawkCalculations::~TomahawkCalculations()
{
}

bool TomahawkCalculations::loadGroups(const std::string& file){
	if(file.size() == 0){
		std::cerr << "No file set" << std::endl;
		return false;
	}

	std::ifstream stream(file);
	if(!stream.good()){
		std::cerr << "bad file" << std::endl;
		return false;
	}

	std::string line;
	U32 n_lines = 1;

	//
	if(!getline(stream, line)){
		std::cerr << "failed to get first line" << std::endl;
		return false;
	}
	std::istringstream ss(line);
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

	std::cerr << "header lines" << std::endl;

	while(getline(stream, line)){
		// Empty lines
		if(line.size() == 0)
			break;

		std::istringstream ss(line);

		// Assert correct format
		// Count tabs until out-of-range
		size_t pos = -1;
		U32 inner_tabs = 0;
		while(true){
			pos = line.find('\t', pos + 1);
			if(pos == std::string::npos)
				break;

			++inner_tabs;
		}

		if(tabs != inner_tabs){
			std::cerr << Helpers::timestamp("ERROR") << "Illegal format! Expected " << tabs << " columns! Line: " << n_lines << "..." << std::endl;
			return false;
		}

		// do stuff
		++n_lines;
	}
	// Get sample names mappings
	// this->totempole.sampleHashTable->GetItem()

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
	case 1: return(this->__calculateSFS<BYTE>());
	case 2: return(this->__calculateSFS<U16>());
	case 4: return(this->__calculateSFS<U32>());
	case 8: return(this->__calculateSFS<U64>());
	default: exit(1); break;
	}

	return false;
}

} /* namespace Tomahawk */
