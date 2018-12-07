#include <sys/time.h>
#include <regex>

#include "utility.h"

namespace tomahawk {
namespace utility {

std::string remove_whitespace(std::string& string){
	string.erase(remove_if(string.begin(), string.end(), isspace), string.end());
	return(string);
}

std::string remove_excess_whitespace(const std::string& string){
	return(std::regex_replace(string, std::regex("^ +| +$|( ) +"), std::string("$1")));
}

int IsBigEndian(){
	union {
		uint32_t i;
		uint8_t  c[4];
	} val = {0x01020304};

	return val.c[0] == 1;
}

std::vector<std::string> &split(std::string& s, char delim, std::vector<std::string>& elems, const bool keepEmpty) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
    	if(!item.empty() || (item.empty() && keepEmpty))
    		elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(std::string& s, char delim, const bool keepEmpty) {
	std::vector<std::string> elems;
	split(s, delim, elems, keepEmpty);
	return elems;
}

std::vector<std::string> &split(const std::string& s, char delim, std::vector<std::string>& elems, const bool keepEmpty) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
    	if(!item.empty() || (item.empty() && keepEmpty))
    		elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string& s, char delim, const bool keepEmpty) {
	std::vector<std::string> elems;
	split(s, delim, elems, keepEmpty);
	return elems;
}

std::vector<std::string> splitLastOf(const std::string& s, const char delim, const bool includeLast){
	std::size_t found = s.find_last_of(delim);
	std::vector<std::string> ret;
	ret.push_back(s.substr(0, found + includeLast));
	ret.push_back(s.substr(found + 1));
	return(ret);
}


std::string datetime(){
	time_t t = time(0);
	struct timeval  tv;
	struct timezone tz;
	struct tm      	*now = localtime(&t);
	gettimeofday(&tv, &tz);

	char buffer[48];
	sprintf(buffer, "%04u-%02u-%02u %02u:%02u:%02u,%03u",
			now->tm_year + 1900,
			now->tm_mon + 1,
			now->tm_mday,
			now->tm_hour,
			now->tm_min,
			now->tm_sec,
			(uint32_t)tv.tv_usec / 1000);

	return std::string(&buffer[0], 23);
}

std::string timestamp(const std::string type){
	std::stringstream ret;

	ret << "[" << datetime() << "]";
	ret << "[" << type << "] ";

	return(ret.str());
}

std::string timestamp(const std::string type, const std::string type2){
	std::stringstream ret;

	ret << "[" << datetime() << "]";
	ret << "[" << type << "]";
	ret << "[" << type2 << "] ";

	return(ret.str());
}

std::string SecondsToTimestring(const double& value){
	uint32_t internalVal = value;
	std::string retVal;
	const uint32_t hours = internalVal / 3600;
	if(hours > 0) retVal += std::to_string(hours) + "h";
	internalVal %= 3600;
	const uint32_t min = internalVal / 60;
	if(min > 0) retVal += std::to_string(min) + "m";
	internalVal %= 60;
	const uint32_t sec = internalVal;
	retVal += std::to_string(sec) + "s";

	return(retVal);
}

std::string BasePath(const std::string& input){
	size_t found = input.find_last_of("/\\");
	if(found != std::string::npos)
		return(input.substr(0, found));
	else return std::string();
}

std::string BaseName(const std::string& input){
	size_t found = input.find_last_of("/\\");
	if(found == std::string::npos)
		found = -1;

	return(input.substr(found+1, input.size()));
}

std::string ExtensionName(const std::string& input){
	std::string base = BaseName(input);
	size_t foundDot = base.rfind('.');
	if(foundDot == std::string::npos)
		foundDot = base.size() - 1;

	return(base.substr(foundDot+1, base.size()));
}


std::vector<std::string> FilePathBaseExtension(const std::string& input){
	std::vector<std::string> ret;

	const std::string base = BaseName(input);
	ret.push_back(BasePath(input));
	ret.push_back(base);

	size_t foundDot = base.rfind('.');
	if(foundDot == std::string::npos){
		ret.push_back(base);
		ret.push_back(std::string());
	} else {
		ret.push_back(base.substr(0, foundDot));
		ret.push_back(base.substr(foundDot+1, base.size()));
	}

	ret.push_back(ExtensionName(input));
	return(ret);
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

int32_t ConvertCharToInt(const char& input){
	if(input >= '0' && input <= '9') return input - '0';
	else if(input >= 'A' && input <= 'F') return input - 'A' + 10;
	else if(input >= 'a' && input <= 'f') return input - 'a' + 10;
	throw std::invalid_argument("Invalid input string");
}

bool HexToBytes(const std::string& hex, uint8_t* target){
	if(hex.size() % 2 != 0){
		std::cerr << "illegal uneven hex" << std::endl;
		return false;
	}

	uint32_t p = 0;
	for (uint32_t i = 0; i < hex.length(); i += 2, ++p)
		target[p] = ConvertCharToInt(hex[i])*16 + ConvertCharToInt(hex[i+1]);

	return true;
}

}

void SerializeString(const std::string& string, std::ostream& stream){
	size_t size_helper = string.size();
	stream.write((const char*)&size_helper, sizeof(size_t));
	stream.write(string.data(), string.size());
}

void DeserializeString(std::string& string, std::istream& stream){
	size_t l_string;
	stream.read((char*)&l_string, sizeof(size_t));
	string.resize(l_string);
	stream.read(&string[0], l_string);
}

}
