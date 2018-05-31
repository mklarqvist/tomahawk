#ifndef TOMAHAWKIMPORTERFILTERS_H_
#define TOMAHAWKIMPORTERFILTERS_H_

namespace tomahawk{

struct ImporterFilters{
	ImporterFilters() :
		dropUnivariantRef(false),
		dropUnivariantAlt(false),
		flipMajorMinor(false),
		MAF(0),
		HWE_P(0),
		missingness(0.2)
	{}
	~ImporterFilters(){}

	bool dropUnivariantRef;
	bool dropUnivariantAlt;
	bool flipMajorMinor;
	double MAF;
	double HWE_P;
	double missingness;
};

}

#endif /* TOMAHAWKIMPORTERFILTERS_H_ */
