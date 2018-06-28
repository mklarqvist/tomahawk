
#include <tomahawk/tomahawk_calc.h>

namespace tomahawk {

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

bool TomahawkCalc::addRegions(std::vector<std::string>& positions){
	return(this->reader.addRegions(positions));
}

template <typename K, typename V>
bool comparePairs(const std::pair<K,V>& a, const std::pair<K,V>& b){ return a.first < b.first; }


bool TomahawkCalc::Calculate(pair_vector& blocks){
	if((this->parameters_validated == false) && (!this->parameters.Validate()))
		return false;

	this->parameters_validated = true;

	std::sort(blocks.begin(), blocks.end(), comparePairs<U32, U32>);
	if(!this->reader.getBlocks(blocks)){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to get Tomahawk blocks..." << std::endl;
		return false;
	}

	const BYTE bit_width = this->reader.getBitWidth();
	if(bit_width == 1) 	    return(this->Calculate<BYTE>());
	else if(bit_width == 2) return(this->Calculate<U16>());
	else if(bit_width == 4) return(this->Calculate<U32>());
	else if(bit_width == 8) return(this->Calculate<U64>());
	else {
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Impossible bit width..." << std::endl;
		exit(1);
	}

	return false;
}

bool TomahawkCalc::Calculate(std::vector<U32>& blocks){
	if((this->parameters_validated == false) && (!this->parameters.Validate()))
		return false;

	this->parameters_validated = true;

	if(!this->reader.getBlocks(blocks)){
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Failed to get Tomahawk blocks..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << helpers::timestamp("LOG","TOMAHAWK") << "Inflated " << blocks.size() << " blocks..." << std::endl;

	const BYTE bit_width = this->reader.getBitWidth();
	if(bit_width == 1) 	    return(this->Calculate<BYTE>());
	else if(bit_width == 2) return(this->Calculate<U16>());
	else if(bit_width == 4) return(this->Calculate<U32>());
	else if(bit_width == 8) return(this->Calculate<U64>());
	else {
		std::cerr << helpers::timestamp("ERROR", "TOMAHAWK") << "Impossible bit width..." << std::endl;
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

	if(this->parameters.window_mode){
		if(reader.interval_tree_entries != nullptr){
			std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Window mode when slicng is not supported..." << std::endl;
			return false;
		}

		if(!this->balancer.BuildWindow(this->reader, this->parameters.n_threads, this->parameters.n_window_bases)){
			std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Failed to split into blocks..." << std::endl;
			return false;
		}
	} else {
		if(!this->balancer.Build(this->reader, this->parameters.n_threads)){
			std::cerr << helpers::timestamp("ERROR", "BALANCER") << "Failed to split into blocks..." << std::endl;
			return false;
		}
	}

	return(this->Calculate(this->balancer.getLoad()));
}

} /* namespace Tomahawk */
