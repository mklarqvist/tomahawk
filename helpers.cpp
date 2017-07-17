/*
POLYGON - Next generation sequencing quality control and preprocessing
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

#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <sys/time.h>

#include "TypeDefinitions.h"
#include "helpers.h"

namespace Tomahawk{
namespace Helpers{

std::vector<std::string> &split(std::string& s, char delim, std::vector<std::string>& elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
    	if (!item.empty())
    		elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(std::string& s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

std::vector<std::string> &split(const std::string& s, char delim, std::vector<std::string>& elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
    	if (!item.empty())
    		elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string& s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

std::vector<std::string> splitLastOf(const std::string& s, const char delim, const bool includeLast){
	std::size_t found = s.find_last_of(delim);
	std::vector<std::string> ret;
	ret.push_back(s.substr(0, found + includeLast));
	ret.push_back(s.substr(found + 1));
	return(ret);
}

std::string output_separator(const U32 width){
	std::string ret;
	ret.reserve(width);
	for(U32 i = 0; i < width; ++i)
		ret += '-';

	return ret;
}

std::string datetime(){
	time_t t = time(0);
	struct timeval  	tv;
	struct timezone 	tz;
	struct tm      	*now = localtime(&t);
	gettimeofday(&tv, &tz);
	std::string sec = std::to_string(now->tm_sec);
	if(sec.size() == 1) sec = '0' + sec;

	const S32 remainder = tv.tv_usec / 1000;
	char pad = 0, pad2 = 0;
	if(remainder < 10){ pad = '0'; pad2 = '0';}
	else if(remainder < 100){pad = '0';}

	std::stringstream ret;
	ret << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-'
		<<  now->tm_mday << " " << now->tm_hour << ":" << now->tm_min << ":"
		<< sec << "," << pad << pad2 << remainder;

		return ret.str();
}

std::string timestamp(const std::string type){
	std::stringstream ret;
	ret << "\33[2K\r\033[0m";
	if(type == "ERROR") ret << "\033[0;31m";
	else if(type == "WARNING") ret << "\033[38;5;208m";

	ret << "[" << datetime() << "]";
	ret << "[" << type << "] ";

	return(ret.str());
}

std::string timestamp(const std::string type, const std::string type2){
	std::stringstream ret;
	//ret << "\33[2K\r\033[0m";
	//if(type == "ERROR" || type2 == "ERROR") ret << "\033[0;31m";
	//else if(type == "WARNING") ret << "\033[38;5;226m";

	ret << "[" << datetime() << "]";
	ret << "[" << type << "]";
	ret << "[" << type2 << "] ";

	return(ret.str());
}

double diffclock(clock_t clock1, clock_t clock2){
	double diffticks = clock1 - clock2;
	double diffms    = diffticks / ( CLOCKS_PER_SEC / 1000 );
	return diffms;
}

std::string NumberThousandsSeparator(std::string number){
	int insertPosition = number.length() - 3;
	char EndPos = 0;
	if(number[0] == '-'){
		EndPos = 1;
	}

	// Todo: fix NNNN.MMMM

	while (insertPosition > EndPos) {
	    number.insert(insertPosition, ",");
	    insertPosition -= 3;
	}

	if(number[0] == '-'){
		std::string numberTemp = "&ndash;";
		numberTemp += number.substr(1,10000);
		return numberTemp;
	}

	return number;
}

std::string MillisecondsToTimestring(const U32 ms){
	const double seconds = ms / 1000;
	const S32 hours = ((S32)seconds / 60 / 60);
	const S32 minutes = ((S32)seconds / 60) % 60;
	const S32 sec = (S32)seconds % 60;
	const S32 remainder = (ms % 1000);

	char pad = 0, pad2 = 0;
	if(remainder < 10){ pad = '0'; pad2 = '0';}
	else if(remainder < 100){pad = '0';}

	std::stringstream st;

	if(hours > 0) st << hours << 'h';
	if(minutes > 0) st << minutes << 'm';
	st << sec << '.' << pad << pad2 << remainder << 's';

	return st.str();
}

}
}
