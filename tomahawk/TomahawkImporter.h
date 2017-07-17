#ifndef TOMAHAWKIMPORTER_H_
#define TOMAHAWKIMPORTER_H_

#include <string>
#include "../TypeDefinitions.h"
#include "../io/reader.h"
#include "../io/vcf/VCFParser.h"


namespace Tomahawk {

class TomahawkImporter {
public:
	TomahawkImporter(std::string inputFile, std::string outputPrefix) : reader_(inputFile), parser_(this->reader_, outputPrefix){}
	virtual ~TomahawkImporter(){}

	bool Open(void);
	bool Build();

private:
	Tomahawk::reader reader_;
	Tomahawk::VCF::VCFParser parser_;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKIMPORTER_H_ */
