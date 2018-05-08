
#include "TomahawkCalc.h"

namespace Tomahawk {

TomahawkCalc::TomahawkCalc(void) :
	parameters_validated(false)
{}

TomahawkCalc::~TomahawkCalc(){}

bool TomahawkCalc::Open(const std::string input, const std::string output){
	if(!this->reader.open(input)){
		return false;
	}

	this->input_file = input;
	this->output_file = output;

	return true;
}

template <typename K, typename V>
bool comparePairs(const std::pair<K,V>& a, const std::pair<K,V>& b){ return a.first < b.first; }


bool TomahawkCalc::Calculate(pair_vector& blocks){
	if((this->parameters_validated == false) && (!this->parameters.Validate()))
		return false;

	this->parameters_validated = true;

	std::sort(blocks.begin(), blocks.end(), comparePairs<U32, U32>);
	if(!this->reader.getBlocks(blocks)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to get Tomahawk blocks..." << std::endl;
		return false;
	}

	const BYTE bit_width = this->reader.getBitWidth();
	if(bit_width == 1) 	    return(this->Calculate<BYTE>());
	else if(bit_width == 2) return(this->Calculate<U16>());
	else if(bit_width == 4) return(this->Calculate<U32>());
	else if(bit_width == 8) return(this->Calculate<U64>());
	else {
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Impossible bit width..." << std::endl;
		exit(1);
	}

	return false;
}

bool TomahawkCalc::Calculate(std::vector<U32>& blocks){
	if((this->parameters_validated == false) && (!this->parameters.Validate()))
		return false;

	this->parameters_validated = true;

	if(!this->reader.getBlocks(blocks)){
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to get Tomahawk blocks..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG","TOMAHAWK") << "Inflated " << blocks.size() << " blocks..." << std::endl;

	const BYTE bit_width = this->reader.getBitWidth();
	if(bit_width == 1) 	    return(this->Calculate<BYTE>());
	else if(bit_width == 2) return(this->Calculate<U16>());
	else if(bit_width == 4) return(this->Calculate<U32>());
	else if(bit_width == 8) return(this->Calculate<U64>());
	else {
		std::cerr << Helpers::timestamp("ERROR", "TOMAHAWK") << "Impossible bit width..." << std::endl;
		exit(1);
	}

	return false;
}

bool TomahawkCalc::Calculate(){
	if((this->parameters_validated == false) && (!this->parameters.Validate()))
		return false;

	this->parameters_validated = true;
	this->balancer.setSelected(this->parameters.chunk_selected);
	this->balancer.setDesired(this->parameters.n_chunks);

	if(!this->balancer.Build(this->reader, this->parameters.n_threads)){
		std::cerr << Helpers::timestamp("ERROR", "BALANCER") << "Failed to split into blocks..." << std::endl;
		return false;
	}

	return(this->Calculate(this->balancer.getLoad()));
}

} /* namespace Tomahawk */
