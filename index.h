#ifndef INDEX_H_
#define INDEX_H_


#include "tomahawk/TomahawkReader.h"
#include "tomahawk/TotempoleReader.h"
#include "io/TomahawkOutput/TomahawkOutputReader.h"

int index(int argc, char** argv){
	argc -= 2; argv += 2;
	programMessage();
	std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling index..." << std::endl;

	if(argc < 1){
		std::cerr << argc << std::endl;
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Missing parameters" << std::endl;
		return(1);
	}

	std::string inputFile(&argv[0][0]);

	// Parse file suffix
	std::vector<std::string> paths = Tomahawk::Helpers::splitLastOf(inputFile, '/', true);
	std::vector<std::string> files = Tomahawk::Helpers::splitLastOf(paths[1], '.');

	// Todo: if failed to read from file suffix: try to look into file header MAGIC
	if(files[1].size() == 0){
		std::cerr << "could not determine file type from suffix" << std::endl;
		return false;
	}

	std::transform(files[1].begin(), files[1].end(), files[1].begin(), ::tolower);

	if(files[1] == Tomahawk::Constants::OUTPUT_SUFFIX){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR","INDEX") << "Twk files are already indexed..." << std::endl;
	} else if(files[1] == Tomahawk::Constants::OUTPUT_LD_SUFFIX) {
		Tomahawk::IO::TomahawkOutputReader reader;

		if(!reader.index(inputFile)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "INDEX") << "Failed to index file!" << std::endl;
			return 1;
		}
	} else {
		std::cerr << "Unknown file type" << std::endl;
	}

	return 0;
}



#endif /* INDEX_H_ */
