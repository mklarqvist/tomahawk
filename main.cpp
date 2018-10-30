#include <fstream>

#include "tomahawk.h"
std::string tomahawk::LITERAL_COMMAND_LINE;
std::string tomahawk::INTERPRETED_COMMAND;

#include "importer.h"
#include "ld.h"
#include "two_reader.h"

#include "import.h"
#include "calc.h"
#include "view.h"

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

	if(strncmp(&argv[1][0], "import", 5) == 0){
		return(import(argc, argv));

	}

	else if(strncmp(&argv[1][0], "calc", 4) == 0){
		return(calc(argc, argv));

	}

	else if(strncmp(&argv[1][0], "view", 4) == 0){
		return(view(argc, argv));

	}
	/*
	else if(strncmp(&argv[1][0], "sort", 4) == 0){
		return(sort(argc, argv));

	} else if(strncmp(&argv[1][0], "concat", 6) == 0){
		return(concat(argc, argv));

	} else if(strncmp(&argv[1][0], "stats", 5) == 0){
		return(stats(argc, argv));

	} else if(strncmp(&argv[1][0], "aggregate", 8) == 0){
		return(aggregate(argc, argv));
	}
	*/
	else if(strncmp(&argv[1][0], "--version", 9) == 0 || strncmp(&argv[1][0], "version", 7) == 0){
		tomahawk::ProgramMessage(false);
		return(0);

	} else if(strncmp(&argv[1][0], "--help", 6) == 0 || strncmp(&argv[1][0], "help", 4) == 0){
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
