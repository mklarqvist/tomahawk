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
