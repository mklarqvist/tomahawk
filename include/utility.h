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
#ifndef TOMAHAWK_UTILS_H_
#define TOMAHAWK_UTILS_H_

#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <iostream>
#include <cstring>
#include <regex>

namespace tomahawk {
namespace utility {

int IsBigEndian(void);

std::vector<std::string> &split(std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(std::string &s, char delim, const bool keepEmpty = false);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim, const bool keepEmpty = false);
std::vector<std::string> splitLastOf(const std::string& s, const char delim, const bool includeLast = false);

std::string remove_whitespace(std::string& string);
std::string remove_excess_whitespace(const std::string& string);

std::string timestamp(const std::string type);
std::string timestamp(const std::string type, const std::string type2);
std::string datetime();
std::string NumberThousandsSeparator(std::string number);

std::string BasePath(const std::string& input);
std::string BaseName(const std::string& input);
std::string ExtensionName(const std::string& input);
std::vector<std::string> FilePathBaseExtension(const std::string& input);

template <class T>
std::string ToPrettyString(const T& data){
	return utility::NumberThousandsSeparator(std::to_string(data));
}

template <class T>
std::string ToPrettyString(const std::vector<T>& data){
	std::string ret;
	for(uint32_t i = 0; i < data.size() - 1; ++i){
		ret += utility::NumberThousandsSeparator(std::to_string(data[i]));
		ret += ", ";
	}
	ret += std::to_string(data[data.size()-1]);
	return ret;
}

std::string SecondsToTimestring(const double& value);

int32_t ConvertCharToInt(const char& input);
bool HexToBytes(const std::string& hex, uint8_t* target);

template <class T>
std::string ToPrettyDiskString(const T value){
	if(value == 0) return("0 b");
	if(value > 1E9){
		return(std::to_string((double)value/1e9) + " Gb");
	} else if(value > 1E6){
		return(std::to_string((double)value/1e6) + " Mb");
	} else if(value > 1E3){
		return(std::to_string((double)value/1e3) + " Kb");
	} else
		return(std::to_string(value) + " b");
}

}

void SerializeString(const std::string& string, std::ostream& stream);
void DeserializeString(std::string& string, std::istream& stream);

template <class T>
void SerializePrimitive(const T& value, std::ostream& stream){
	stream.write((const char*)&value, sizeof(T));
}

template <class T>
void DeserializePrimitive(T& value, std::istream& stream){
	stream.read((char*)&value, sizeof(T));
}

}

#endif
