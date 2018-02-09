#ifndef TOMAHAWKIMPORTERFILTERS_H_
#define TOMAHAWKIMPORTERFILTERS_H_

namespace Tomahawk{

struct ImporterFilters{
	ImporterFilters() : MAF(0), HWE_P(0), missingness(0.2){}
	~ImporterFilters(){}

	double MAF;
	double HWE_P;
	double missingness;
};

}

#endif /* TOMAHAWKIMPORTERFILTERS_H_ */
