#include "TomahawkCalculations.h"

namespace Tomahawk {

TomahawkCalculations::TomahawkCalculations()
{
}

TomahawkCalculations::~TomahawkCalculations()
{
}

bool TomahawkCalculations::calculateTajimaD(const U32 bin_size){
	if(this->Occ.size() == 0){
		switch(this->bit_width){
		case 1: return(this->__calculateTajimaD<BYTE>(bin_size));
		case 2: return(this->__calculateTajimaD<U16>(bin_size));
		case 4: return(this->__calculateTajimaD<U32>(bin_size));
		case 8: return(this->__calculateTajimaD<U64>(bin_size));
		default: exit(1); break;
		}
	} else {
		switch(this->bit_width){
		case 1: return(this->__calculateTajimaDGrouped<BYTE>(bin_size));
		case 2: return(this->__calculateTajimaDGrouped<U16>(bin_size));
		case 4: return(this->__calculateTajimaDGrouped<U32>(bin_size));
		case 8: return(this->__calculateTajimaDGrouped<U64>(bin_size));
		default: exit(1); break;
		}
	}

	return false;
}

bool TomahawkCalculations::calculateSiteStats(void){
	if(this->Occ.size() == 0){
		switch(this->bit_width){
		case 1: return(this->__calculateSiteStats<BYTE>());
		case 2: return(this->__calculateSiteStats<U16>());
		case 4: return(this->__calculateSiteStats<U32>());
		case 8: return(this->__calculateSiteStats<U64>());
		default: exit(1); break;
		}
	} else {
		switch(this->bit_width){
		case 1: return(this->__calculateSiteStatsGrouped<BYTE>());
		case 2: return(this->__calculateSiteStatsGrouped<U16>());
		case 4: return(this->__calculateSiteStatsGrouped<U32>());
		case 8: return(this->__calculateSiteStatsGrouped<U64>());
		default: exit(1); break;
		}
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
