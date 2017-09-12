#ifndef TOMAHAWK_UTILITY_H_
#define TOMAHAWK_UTILITY_H_

#include <iostream>
#include "support/MagicConstants.h"

// Declare extern
char Tomahawk::Constants::PROGRAM_VERSION_FRONT[5048];
char Tomahawk::Constants::PROGRAM_VERSION_BACK[5048];
std::string Tomahawk::Constants::LITERAL_COMMAND_LINE;
std::string Tomahawk::Constants::INTERPRETED_COMMAND;

void programMessage(void){
	std::cerr << "Program: " << Tomahawk::Constants::PROGRAM_NAME << " " << Tomahawk::Constants::PROGRAM_VERSION_FRONT << std::endl;
	std::cerr << "Contact: Marcus D. R. Klarqvist <mk21@sanger.ac.uk>" << std::endl;
	std::cerr << "Documentation: https://github.com/mklarqvist/Tomahawk" << std::endl;
	std::cerr << "License: MIT" << std::endl;
	std::cerr << "----------" << std::endl;
}

void programHelp(void){
	std::cerr << "Usage: " << Tomahawk::Constants::PROGRAM_NAME << " [command] [-<commands> [options]*]+" << std::endl;
	std::cerr << "Commands: import, view, calc, sort, index, stats, concat" << std::endl;
}

void programHelpDetailed(void){
	programHelp();
	std::cerr <<
    "\n"
	"concat       concatenate TWO files from the same set of samples\n"
	"calc         calculate linkage disequilibrium (TWO/TOI format)\n"
	"import       import VCF/BCF to TWK/TWI\n"
	"stat         basic statistics of TWK/TWO\n"
	"sort         sort TWO file\n"
    "merge        merge TWO files from non-overlapping sample sets\n"
    "view         TWK->VCF conversion, TWO/TWK view, TWK/TWO subset and filter" << std::endl;
}

#endif /* TOMAHAWK_UTILITY_H_ */
