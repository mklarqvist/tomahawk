#ifndef LIB_HEADER_INTERNAL_H_
#define LIB_HEADER_INTERNAL_H_

// htslib dependencies.
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"

#include "header.h"

namespace tomahawk {

// Returns the hrec that contains information or nullptr if none does.
const bcf_hrec_t* GetPopulatedHrec(const bcf_idpair_t& idPair);

class VcfHeaderInternal : public VcfHeader {
public:
	// Adds Contig information from the idPair to the ContigInfo object.
	void AddContigInfo(const bcf_idpair_t& idPair);

	// Adds FILTER information from the bcf_hrec_t to the VcfFilterInfo object.
	void AddFilterInfo(const bcf_hrec_t* hrec);

	// Adds INFO information from the bcf_hrec_t to the VcfInfo object.
	void AddInfo(const bcf_hrec_t* hrec);

	// Adds FORMAT information from the bcf_hrec_t to the VcfFormatInfo object.
	void AddFormatInfo(const bcf_hrec_t* hrec);

	// Adds structured information from the bcf_hrec_t to the VcfStructuredExtra.
	void AddStructuredExtra(const bcf_hrec_t* hrec);

	// Adds unstructured information from the bcf_hrec_t to the VcfExtra object.
	void AddExtra(const bcf_hrec_t* hrec);

	// Adds sample information from a string.
	void AddSample(const std::string& sample_name);

};

}



#endif /* LIB_HEADER_INTERNAL_H_ */
