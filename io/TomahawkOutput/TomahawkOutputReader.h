#ifndef TOMAHAWKOUTPUTREADER_H_
#define TOMAHAWKOUTPUTREADER_H_

#include <vector>
#include <fstream>
#include <string>
#include <stddef.h>
#include <regex>

#include "../../TypeDefinitions.h"
#include "TomahawkOutputEntry.h"
#include "../../support/PackedEntryReader.h"
#include "TomahawkOutputFilterController.h"


namespace Tomahawk {
namespace IO {
//
//
class TomahawkOutputReader {
	typedef TomahawkOutputEntry entry_type;
	typedef TomahawkOutputFilterController filter_type;

public:
	TomahawkOutputReader();
	~TomahawkOutputReader(){ }

	bool view(const std::string& filename);
	bool index(const std::string& filename);
	bool summary(const std::string& input);

	inline bool matchPositionalString(const std::string& param) const{
		return(std::regex_match(param, std::regex(
				"^"
				"([a-zA-Z0-9_\\.\\-]+[\\$]{0,1}){0,1}"
				"([:]{1}[0-9]+"
				"(([eE]{1}[0-9]){1,}){0,1}"
				"([-]{1}[0-9]+"
				"(([eE]{1}[0-9]){1,}){0,1})?"
				")?"
				"$"
		)));
	}
	inline bool parsePositionalString(std::string& param){
		std::size_t found = param.find(',');
		if(found != std::string::npos){
			std::cerr << "found: " << found << '/' << param.size() << std::endl;
			std::vector<std::string> ret = Tomahawk::Helpers::split(param, ',');
			if(ret.size() != 2){
				std::cerr << "illegal format" << std::endl;
				return 1;
			}
			if(!this->matchPositionalString(ret[0])){
				std::cerr << "failed 1 " << std::endl;
				return false;
			}
			if(!this->matchPositionalString(ret[1])){
				std::cerr << "failed 2 " << std::endl;
				return false;
			}
			return true;
		}

		return(this->matchPositionalString(param));
	}

	// Read entire file into memory
	filter_type& getFilter(void){ return this->filter; }

private:
	bool __viewOnly(void);
	bool __viewFilter(void);

public:
	PackedEntryReader<entry_type, sizeof(entry_type)> reader;
	filter_type filter;
};

}
} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTREADER_H_ */
