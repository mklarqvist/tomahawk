/*
Tomahawk - Next generation sequencing quality control and preprocessing
Copyright (C) 2015-2016, Marcus D. R. Klarqvist.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#ifndef HELPERS_H_
#define HELPERS_H_

#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <atomic>
#include <iostream>
#include <cstring>

#include "TypeDefinitions.h"

namespace Tomahawk{
namespace Helpers{

std::vector<std::string> &split(std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(std::string &s, char delim);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<std::string> splitLastOf(const std::string& s, const char delim, const bool includeLast = false);

std::string timestamp(const std::string type);
std::string timestamp(const std::string type, const std::string type2);
std::string output_separator(const U32 width = 75);
std::string datetime();
double diffclock(clock_t clock1, clock_t clock2);
std::string NumberThousandsSeparator(std::string number);

template <class T>
inline T roundUp(T numToRound, int multiple)
{
	if(multiple == 0) return numToRound;

	int remainder = numToRound % multiple;
	if (remainder == 0) return numToRound;
	return numToRound + multiple - remainder;
}

template <class T>
std::string ToString(const T& data){
	return std::to_string(data);
}

template <class T>
std::string ToPrettyString(const T& data){
	return Tomahawk::Helpers::NumberThousandsSeparator(std::to_string(data));
}

template <class T>
std::string ToPrettyString(const std::vector<T>& data){
	std::string ret;
	for(U32 i = 0; i < data.size() - 1; ++i){
		ret += Tomahawk::Helpers::NumberThousandsSeparator(std::to_string(data[i]));
		ret += ", ";
	}
	ret += std::to_string(data[data.size()-1]);
	return ret;
}

bool matchPositionalStringTWO(const std::string& param);
bool parsePositionalStringTWO(std::string& param);


}
}
#endif /* HELPERS_H_ */
