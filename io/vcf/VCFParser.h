#ifndef VCFPARSER_H_
#define VCFPARSER_H_

#include "VCFHeader.h"
#include "../../tomahawk/TomahawkImportWriter.h"

namespace Tomahawk {
namespace VCF{

class VCFParser {
	typedef VCFParser self_type;
	typedef Tomahawk::reader reader_type;
	typedef VCFHeader header_type;
	typedef TomahawkImportWriter writer_type;
	typedef VCFLine line_type;

public:
	VCFParser(std::string inputFile, std::string outputPrefix);
	~VCFParser();
	bool Build();

private:
	std::string outputPrefix;
	reader_type reader_;
	header_type header_;
	writer_type writer_;
};


} /* namespace VCF */
} /* namespace Tomahawk */

#endif /* VCFPARSER_H_ */
