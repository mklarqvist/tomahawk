#ifndef VCFPARSER_H_
#define VCFPARSER_H_

#include "VCFHeader.h"
#include "../../tomahawk/TomahawkImportWriter.h"

namespace Tomahawk {
namespace VCF{

class VariantParser{
public:
	VariantParser();
	virtual ~VariantParser();

	virtual bool getLine(void) =0;
	virtual bool runLengthEncode() =0;
};

class VCFParser {
	typedef VCFParser selfType;
	typedef Tomahawk::reader& readerType;

public:
	VCFParser(readerType reader, const std::string outputPrefix);
	~VCFParser();
	bool Build();

protected:
	bool GetHeaderLines(void);
	bool ValidateVCF(void);
	bool SampleLine(void);

private:
	std::string outputPrefix;
	readerType reader_;
	VCFHeader header_;
	TomahawkImportWriter writer_;
};


} /* namespace VCF */
} /* namespace Tomahawk */

#endif /* VCFPARSER_H_ */
