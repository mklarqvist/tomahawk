#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <regex>

#include "TypeDefinitions.h"
#include "helpers.h"
#include "simd_definitions.h"
#include "MagicConstants.h"

namespace Tomahawk{
namespace Helpers{

int isBigEndian(){
	union {
		uint32_t i;
		uint8_t  c[4];
	} val = {0x01020304};

	return val.c[0] == 1;
}

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

std::string program_string(const bool back){
	return(std::string(Constants::LITERAL_COMMAND_LINE)
			+ "; VERSION=" + std::string(back ? Tomahawk::Constants::PROGRAM_VERSION_BACK : Tomahawk::Constants::PROGRAM_VERSION_FRONT)
			+ "; Date=" + Tomahawk::Helpers::datetime() + "; SIMD=" + SIMD_MAPPING[SIMD_VERSION]);
}

std::string datetime(){
	time_t t = time(0);
	struct timeval  tv;
	struct timezone tz;
	struct tm      	*now = localtime(&t);
	gettimeofday(&tv, &tz);

	char buffer[23];
	sprintf(buffer, "%04u-%02u-%02u %02u:%02u:%02u,%03u",
			now->tm_year + 1900,
			now->tm_mon + 1,
			now->tm_mday,
			now->tm_hour,
			now->tm_min,
			now->tm_sec,
			(U32)tv.tv_usec / 1000);

	return std::string(&buffer[0], 23);
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

	ret << "\33[2K\r\033[0m";
	if(type == "ERROR") ret << "\033[0;31m";
	else if(type == "WARNING") ret << "\033[38;5;208m";

	ret << "[" << datetime() << "]";
	ret << "[" << type << "]";
	ret << "[" << type2 << "] ";

	return(ret.str());
}

std::string NumberThousandsSeparator(std::string number){
	int insertPosition = number.length() - 3;
	char EndPos = 0;
	if(number[0] == '-')
		EndPos = 1;

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

bool matchPositionalStringTWO(const std::string& param){
	return(std::regex_match(param, std::regex(
			"^"
			"([a-zA-Z0-9_\\.\\-\\|]+[\\$]{0,1})+"
			"([:]{1}"
			"([0-9]{1,}([\\.]{1}[0-9]{1,}){0,1}([eE]{1}[0-9]{1,}){0,1}){1}"
			"([-]{1}"
			"([0-9]{1,}([\\.]{1}[0-9]{1,}){0,1}([eE]{1}[0-9]{1,}){0,1}){1})?"
			"){0,1}"
			"$"
	)));
}

bool parsePositionalStringTWO(const std::string& param){
	std::size_t found = param.find(',');
	if(found != std::string::npos){
		std::vector<std::string> ret = Tomahawk::Helpers::split(param, ',');
		if(ret.size() != 2){
			std::cerr << "illegal format" << std::endl;
			return 1;
		}
		if(!matchPositionalStringTWO(ret[0])){
			std::cerr << "failed 1 " << std::endl;
			return false;
		}
		if(!matchPositionalStringTWO(ret[1])){
			std::cerr << "failed 2 " << std::endl;
			return false;
		}
		return true;
	}
	return(matchPositionalStringTWO(param));
}

}
}
