#include <fstream>

#include "tomahawk.h"
std::string tomahawk::LITERAL_COMMAND_LINE;
std::string tomahawk::INTERPRETED_COMMAND;

#include "import.h"
#include "calc.h"
#include "view.h"
#include "concat.h"
#include "sort.h"
#include "aggregate.h"
#include "stats.h"
#include "haplotype.h"
#include "relationship.h"
#include "decay.h"
#include "scalc.h"

int main(int argc, char** argv){
	if(tomahawk::utility::IsBigEndian()){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Tomahawk does not support big endian systems..." << std::endl;
		return(1);
	}

	if(argc == 1){
		tomahawk::ProgramMessage();
		tomahawk::ProgramHelpDetailed();
		return(1);
	}

	// Literal string input line
	tomahawk::LITERAL_COMMAND_LINE = tomahawk::TOMAHAWK_PROGRAM_NAME;
	for(int i = 1; i < argc; ++i)
		tomahawk::LITERAL_COMMAND_LINE += " " + std::string(&argv[i][0]);

	if(strcmp(&argv[1][0], "import") == 0){
		return(import(argc, argv));

	}

	else if(strcmp(&argv[1][0], "calc") == 0){
		return(calc(argc, argv));

	} else if(strcmp(&argv[1][0], "calc-single") == 0 || strcmp(&argv[1][0], "scalc") == 0){
		return(scalc(argc, argv));
	}

	else if(strcmp(&argv[1][0], "view") == 0){
		return(view(argc, argv));

	}

	else if(strncmp(&argv[1][0], "sort", 4) == 0){
		return(sort(argc, argv));

	}

	else if(strncmp(&argv[1][0], "concat", 6) == 0){
		return(concat(argc, argv));

	}
	else if(strncmp(&argv[1][0], "stats", 5) == 0){
		return(stats(argc, argv));

	}
	else if(strncmp(&argv[1][0], "aggregate", 8) == 0){
		return(aggregate(argc, argv));
	}
	else if(strncmp(&argv[1][0], "haplotype", 9) == 0){
		return(haplotype(argc, argv));
	}
	else if(strncmp(&argv[1][0], "relationship", 9) == 0){
		return(relationship(argc, argv));
	}
	else if(strncmp(&argv[1][0], "decay", 5) == 0){
		return(decay(argc, argv));
	}
	else if(strcmp(&argv[1][0], "--version") == 0 || strcmp(&argv[1][0], "version") == 0){
		tomahawk::ProgramMessage(false);
		return(0);

	} else if(strcmp(&argv[1][0], "--help") == 0 || strcmp(&argv[1][0], "help") == 0){
		tomahawk::ProgramMessage();
		tomahawk::ProgramHelpDetailed();
		return(0);

	} else {
		tomahawk::ProgramMessage();
		tomahawk::ProgramHelpDetailed();
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Illegal command" << std::endl;
		return(1);
	}
	return(1);
}
