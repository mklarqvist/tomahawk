#ifndef TOMAHAWK_UTILITY_H_
#define TOMAHAWK_UTILITY_H_

#include <iostream>
#include "support/MagicConstants.h"

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

#endif /* TOMAHAWK_UTILITY_H_ */
