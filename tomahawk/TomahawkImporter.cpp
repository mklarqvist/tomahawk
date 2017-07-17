#include <cstring>

#include "TomahawkImporter.h"

namespace Tomahawk {

bool TomahawkImporter::Open(void){
	if(!this->reader_.open())
		return false;

	return true;
}

bool TomahawkImporter::Build(){
	if(!this->Open())
		return false;

	if(!this->parser_.Build())
		return false;

	return true;
}


} /* namespace Tomahawk */
