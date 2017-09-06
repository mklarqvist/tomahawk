/*
Copyright (C) 2016-2017 Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk21@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/
#include "utility.h"
#include "tomahawk/TomahawkReader.h"
#include "totempole/TotempoleReader.h"
#include "tomahawk/TomahawkOutput/TomahawkOutputReader.h"
#include "algorithm/sort/TomahawkOutputSort.h"

int sort(int argc, char** argv){
	argc -= 2; argv += 2;
	programMessage();
	std::cerr << Tomahawk::Helpers::timestamp("LOG") << "Calling sort..." << std::endl;

	if(argc < 2){
		std::cerr << argc << std::endl;
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Missing parameters" << std::endl;
		return(1);
	}

	std::string inputFile(&argv[0][0]);
	std::string outputFile(&argv[1][0]);

	// is merge?
	bool merge = argv[2];
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
		Tomahawk::Algorithm::Output::TomahawkOutputSorter reader;

		if(!merge){
			if(!reader.sort(inputFile, outputFile, 1e9)){
				std::cerr << Tomahawk::Helpers::timestamp("ERROR", "SORT") << "Failed to sort file!" << std::endl;
				return 1;
			}
		} else {
			std::cerr << "is merge" << std::endl;
			if(!reader.sortMerge(inputFile)){
				std::cerr << "failed merge" << std::endl;
				return 1;
			}

		}
	} else {
		std::cerr << "Unknown file type" << std::endl;
	}

	return 0;
}
