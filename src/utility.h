#ifndef TOMAHAWK_UTILITY_H_
#define TOMAHAWK_UTILITY_H_

#include <iostream>
#include "support/MagicConstants.h"

// Declare extern
std::string Tomahawk::Constants::LITERAL_COMMAND_LINE;
std::string Tomahawk::Constants::INTERPRETED_COMMAND;

void programMessage(const bool separator = true){
	std::cerr << "Program: " << Tomahawk::Constants::PROGRAM_NAME << " " << VERSION << std::endl;
	std::cerr << "Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>" << std::endl;
	std::cerr << "Documentation: https://github.com/mklarqvist/Tomahawk" << std::endl;
	std::cerr << "License: MIT" << std::endl;
	if(separator) std::cerr << "----------" << std::endl;
}

void programHelp(void){
	std::cerr << "Usage: " << Tomahawk::Constants::PROGRAM_NAME << " [--version] [--help] <commands> <argument>" << std::endl;
	std::cerr << "Commands: import, view, calc, sort, index, stats, concat, tajima, fst, abba, sfs" << std::endl;
}

void programHelpDetailed(void){
	programHelp();
	std::cerr <<
    "\n"
	"concat       concatenate TWO files from the same set of samples\n"
	"calc         calculate linkage disequilibrium (TWO/TOI format)\n"
	"import       import VCF/BCF to TWK/TWI\n"
	"stats        basic statistics of TWK/TWO\n"
	"sort         sort TWO file\n"
    "view         TWK->VCF conversion, TWO/TWK view, TWK/TWO subset and filter\n"
	"tajima       calculate Tajima's D, nucleotide diversity, and mean AF\n"
	"sfs          calculate the site frequency spectrum (SFS); a.k.a AFS" << std::endl;
}

#endif /* TOMAHAWK_UTILITY_H_ */
