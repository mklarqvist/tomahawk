#ifndef VCFPARSER_H_
#define VCFPARSER_H_

#include "VCFHeader.h"
#include "../../tomahawk/TomahawkImportWriter.h"

namespace Tomahawk {
namespace VCF{

class VCFParser {
	typedef VCFParser selfType;
	typedef Tomahawk::reader& readerType;
	typedef VCFHeader header_type;
	typedef TomahawkImportWriter writer_type;

public:
	VCFParser(readerType reader, const std::string outputPrefix);
	~VCFParser();
	bool Build();

private:
	std::string outputPrefix;
	readerType reader_;
	VCFHeader header_;
	writer_type writer_;
};


} /* namespace VCF */
} /* namespace Tomahawk */

#endif /* VCFPARSER_H_ */
