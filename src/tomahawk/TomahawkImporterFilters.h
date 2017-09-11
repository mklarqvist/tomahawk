#ifndef TOMAHAWKIMPORTERFILTERS_H_
#define TOMAHAWKIMPORTERFILTERS_H_

namespace Tomahawk{

struct TomahawkImporterFilters{
	TomahawkImporterFilters() : MAF(0), HWE_P(0), missingness(0.2){}
	~TomahawkImporterFilters(){}

	double MAF;
	double HWE_P;
	double missingness;
};

}

#endif /* TOMAHAWKIMPORTERFILTERS_H_ */
