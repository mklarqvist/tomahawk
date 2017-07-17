#ifndef SORT_H_
#define SORT_H_

#include "tomahawk/TomahawkReader.h"
#include "tomahawk/TotempoleReader.h"
#include "io/TomahawkOutput/TomahawkOutputReader.h"
#include "algorithm/sort/TomahawkOutputSort.h"

int sort(int argc, char** argv){
	argc -= 2; argv += 2;
	programMessage();
	std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling sort..." << std::endl;

	if(argc < 1){
		std::cerr << argc << std::endl;
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Missing parameters" << std::endl;
		return(1);
	}

	std::string inputFile(&argv[0][0]);

	// is merge?
	bool merge = argv[1];
	std::cerr << merge << std::endl;

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
		std::cerr << Tomahawk::Helpers::timestamp("ERROR","SORT") << "Twk files are guaranteed sorted..." << std::endl;
	} else if(files[1] == Tomahawk::Constants::OUTPUT_LD_SUFFIX) {
		Tomahawk::Algorithm::Output::TomahawkOutputSorter<Tomahawk::IO::TomahawkOutputEntry> reader;

		if(!merge){
			if(!reader.sort(inputFile)){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR", "SORT") << "Failed to sort file!" << std::endl;
				return 1;
			}
		} else {
			std::cerr << "is merge" << std::endl;
			if(!reader.kwayMerge(inputFile)){
				std::cerr << "failed merge" << std::endl;
				return 1;
			}

		}
	} else {
		std::cerr << "Unknown file type" << std::endl;
	}

	return 0;
}

#endif /* SORT_H_ */
