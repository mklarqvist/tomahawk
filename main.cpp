#include <fstream>

#include "importer.h"

int main(void){

	tomahawk::twk_variant_importer importer;
	return(importer.Import() != 0);
}
