/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

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
==============================================================================*/
#ifndef TOMAHAWK_H_
#define TOMAHAWK_H_

#include <iostream>
#include <string>
#include <regex>

extern int SILENT;

namespace tomahawk {

extern std::string LITERAL_COMMAND_LINE;
extern std::string INTERPRETED_COMMAND;

/*------   Version   ------*/
const int32_t TOMAHAWK_VERSION_MAJOR = 0;
const int32_t TOMAHAWK_VERSION_MINOR = 7;
const int32_t TOMAHAWK_VERSION_PATCH = 0;
const int32_t TOMAHAWK_VERSION_NUMBER  = (TOMAHAWK_VERSION_MAJOR *100*100 + TOMAHAWK_VERSION_MINOR *100 + TOMAHAWK_VERSION_PATCH);
const std::string TOMAHAWK_LIB_VERSION = std::to_string(TOMAHAWK_VERSION_MAJOR) + '.' + std::to_string(TOMAHAWK_VERSION_MINOR) + '.' + std::to_string(TOMAHAWK_VERSION_PATCH);

/*------   Basics   ------*/
const std::string TOMAHAWK_PROGRAM_NAME  = "tomahawk";
const std::string TOMAHAWK_OUTPUT_SUFFIX = "twk";
const std::string TOMAHAWK_MAGIC_HEADER  = "TOMAHAWK\1";
const uint32_t    TOMAHAWK_MAGIC_HEADER_LENGTH = 9;
const std::string TOMAHAWK_LD_SUFFIX = "two";
const std::string TOMAHAWK_LD_MAGIC_HEADER  = "TWO\1";
const uint32_t    TOMAHAWK_LD_MAGIC_HEADER_LENGTH = 4;

/*------   Regular expression patterns  ------*/
const std::regex TWK_REGEX_CANONICAL_BASES = std::regex("^[ATGC]{1}$");
const std::regex TWK_REGEX_CONTIG_ONLY     = std::regex("^[A-Za-z0-9\\-_]+$");
const std::regex TWK_REGEX_CONTIG_POSITION = std::regex("^[A-Za-z0-9\\-_]+\\:[0-9]+([\\.]{1}[0-9]+){0,1}([eE]{1}[0-9]{1})?$");
const std::regex TWK_REGEX_CONTIG_RANGE    = std::regex("^[A-Za-z0-9\\-_]+\\:[0-9]+([\\.]{1}[0-9]+){0,1}([eE]{1}[0-9]{1})?\\-[0-9]+([\\.]{1}[0-9]+){0,1}([eE]{1}[0-9]{1})?$");
const std::regex TWK_REGEX_FLOATING        = std::regex("^[-+]?[0-9]*\\.?[0-9]+$");
const std::regex TWK_REGEX_FLOATING_EXP    = std::regex("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$");
const std::regex TWK_REGEX_NUMBER          = std::regex("^[0-9]+$");
const std::regex TWK_REGEX_NUMBER_EXP      = std::regex("^[0-9]+([eE][-+]?[0-9]+)?$");

/*------   EOF markers   ------*/
const std::string TOMAHAWK_FILE_EOF = "a4f54f39f5e251a6993796f48164ccf554f1b680c2ebbb13be301f3ff76f82cf";
const uint32_t    TOMAHAWK_FILE_EOF_LENGTH = 32;
const uint64_t    TOMAHAWK_INDEX_START_MARKER = 1954702206512158641;

}

#endif
