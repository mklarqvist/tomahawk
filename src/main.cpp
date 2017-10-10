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
#include <iostream>
#include <getopt.h>

#include "utility.h"
#include "calc.h"
#include "import.h"
#include "view.h"
#include "sort.h"
#include "index.h"
#include "concat.h"
#include "stats.h"
#include "tajima.h"
#include "sfs.h"

int main(int argc, char** argv){
	if(Tomahawk::Helpers::isBigEndian()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Tomahawk does not support big endian systems..." << std::endl;
		return(1);
	}

	if(argc == 1){
		programMessage();
		programHelpDetailed();
		return(1);
	}

	// Literal string input line
	Tomahawk::Constants::LITERAL_COMMAND_LINE = Tomahawk::Constants::PROGRAM_NAME;
	for(U32 i = 1; i < argc; ++i)
		Tomahawk::Constants::LITERAL_COMMAND_LINE += " " + std::string(&argv[i][0]);

	if(strncmp(&argv[1][0], "import", 5) == 0){
		return(import(argc, argv));

	} else if(strncmp(&argv[1][0], "calc", 4) == 0){
		return(calc(argc, argv));

	} else if(strncmp(&argv[1][0], "view", 4) == 0){
		return(view(argc, argv));

	} else if(strncmp(&argv[1][0], "sort", 4) == 0){
		return(sort(argc, argv));

	} else if(strncmp(&argv[1][0], "index", 5) == 0){
		//return(index(argc, argv));
		std::cerr << "Not implemented" << std::endl;
		return(1);

	} else if(strncmp(&argv[1][0], "concat", 6) == 0){
		return(concat(argc, argv));

	} else if(strncmp(&argv[1][0], "stats", 5) == 0){
		return(stats(argc, argv));

	} else if(strncmp(&argv[1][0], "tajima", 6) == 0){
		return(tajida(argc, argv));

	} else if(strncmp(&argv[1][0], "sfs", 3) == 0){
		return(sfs(argc, argv));

	} else if(strncmp(&argv[1][0], "--version", 9) == 0 || strncmp(&argv[1][0], "version", 7) == 0){
		programMessage(false);
		return(0);

	} else if(strncmp(&argv[1][0], "--help", 6) == 0 || strncmp(&argv[1][0], "help", 4) == 0){
		programMessage();
		programHelpDetailed();
		return(0);

	} else {
		programMessage();
		programHelpDetailed();
		std::cerr << Tomahawk::Helpers::timestamp("ERROR") << "Illegal command" << std::endl;
		return(1);
	}
	return(1);
}
