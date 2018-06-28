#ifndef TOMAHAWK_UTILITY_H_
#define TOMAHAWK_UTILITY_H_

#include <iostream>

#include "support/magic_constants.h"

// Declare extern
std::string tomahawk::constants::LITERAL_COMMAND_LINE;
std::string tomahawk::constants::INTERPRETED_COMMAND;

void programMessage(const bool separator = true){
	std::cerr << "Program: " << tomahawk::constants::PROGRAM_NAME << " " << VERSION << std::endl;
	std::cerr << "Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>" << std::endl;
	std::cerr << "Documentation: https://github.com/mklarqvist/tomahawk" << std::endl;
	std::cerr << "License: MIT" << std::endl;
	if(separator) std::cerr << "----------" << std::endl;
}

void programHelp(void){
	std::cerr << "Usage: " << tomahawk::constants::PROGRAM_NAME << " [--version] [--help] <commands> <argument>" << std::endl;
	std::cerr << "Commands: import, view, calc, sort, concat, aggregate" << std::endl;
}

void programHelpDetailed(void){
	programHelp();
	std::cerr <<
    "\n"
	"calc         calculate linkage disequilibrium\n"
	"concat       concatenate TWO files from the same set of samples\n"
	"import       import VCF/BCF to TWK\n"
	"sort         sort TWO file\n"
    "view         TWK->VCF conversion, TWO/TWK view, TWK/TWO subset and filter\n"
    "aggregate    data rasterization framework for TWO files\n"<< std::endl;
}

#endif /* TOMAHAWK_UTILITY_H_ */
