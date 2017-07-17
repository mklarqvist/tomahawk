#include "utility.h"
#include "tomahawk/TomahawkImporter.h"

int import(int argc, char** argv){
	argc -= 2; argv += 2;
	programMessage();
	std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling import..." << std::endl;
	if(argc < 2){
		std::cerr << argc << std::endl;
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Missing parameters" << std::endl;
		return(1);
	}
	std::string inputFile(&argv[0][0]);
	std::string output(&argv[1][0]);

	Tomahawk::TomahawkImporter importer(inputFile, output);
	if(!importer.Build()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Failed build!" << std::endl;
		return 1;
	}

	return 0;
}
