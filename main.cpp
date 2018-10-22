#include <fstream>

#include "tomahawk.h"
std::string tomahawk::LITERAL_COMMAND_LINE;
std::string tomahawk::INTERPRETED_COMMAND;


#include "importer.h"

int main(void){

	tomahawk::twk_variant_importer importer;
	return(importer.Import() != 0);
}
