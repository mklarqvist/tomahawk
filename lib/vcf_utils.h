#ifndef TWK_VCFUTILS_H_
#define TWK_VCFUTILS_H_

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

// htslib dependencies.
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"

#include "buffer.h"
#include "support_vcf.h"

namespace tomahawk {
namespace io {

// Returns the hrec that contains information or nullptr if none does.
const bcf_hrec_t* GetPopulatedHrec(const bcf_idpair_t& idPair);

class VcfHeader {
public:
	typedef VcfHeader self_type;
	typedef VcfContig contig_type;
	typedef bcf_hdr_t hts_vcf_header;
	typedef VcfFormat format_type;
	typedef VcfInfo   info_type;
	typedef VcfFilter filter_type;
	typedef VcfStructuredExtra structured_extra_type;
	typedef VcfExtra  extra_type;
	typedef std::unordered_map<std::string, uint32_t> map_type;
	typedef std::unordered_map<uint32_t, uint32_t>    map_reverse_type;

public:
	VcfHeader() = default;
	VcfHeader(const VcfHeader& other);
	~VcfHeader() = default;

	inline size_t GetNumberSamples(void) const{ return(this->samples_.size()); }
	inline size_t GetNumberContigs(void) const{ return(this->contigs_.size()); }

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

	VcfContig* GetContig(const std::string& name);
	VcfContig* GetContig(const int& idx);
	VcfInfo* GetInfo(const std::string& name);
	VcfInfo* GetInfo(const int& idx);
	VcfFormat* GetFormat(const std::string& name);
	VcfFormat* GetFormat(const int& idx);
	VcfFilter* GetFilter(const std::string& name);
	VcfFilter* GetFilter(const int& idx);
	std::string* GetSample(const std::string& name);

	const VcfContig* GetContig(const std::string& name) const;
	const VcfContig* GetContig(const int& idx) const;
	const VcfInfo* GetInfo(const std::string& name) const;
	const VcfInfo* GetInfo(const int& idx) const;
	const VcfFormat* GetFormat(const std::string& name) const;
	const VcfFormat* GetFormat(const int& idx) const;
	const VcfFilter* GetFilter(const std::string& name) const;
	const VcfFilter* GetFilter(const int& idx) const;
	const std::string* GetSample(const std::string& name) const;

	bool BuildReverseMaps(void);
	bool BuildMaps(void);

	/**<
	* Converts this header object into a hts_vcf_header object from the
	* internally stored literal string. This object is required for
	* writing out VCF/BCF files using htslib.
	* @return Returns a bcf_hdr_t pointer.
	*/
	hts_vcf_header* ConvertVcfHeader(void);

	// Append a string to the literal string
	inline void AppendLiteralString(const std::string& literal_addition){ this->literals_ += literal_addition; }

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const self_type& self);
	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, self_type& self);

public:
	// VCF file version string.
	std::string fileformat_string_;
	// Literal string for VcfHeader data. Contains all of the Vcf header data up
	// to the start of the main header line ("#CHROM"...). As such, sample names
	// are not available in this string and needs to be appended before converting
	// back into a htslib vcf header.
	std::string literals_;

	// Vcf header lines parse into:
	// Samples:   Individual sample names.
	// VcfContig: Information relating to the interpretation of a contig. Data
	//            include its name, length in bases, its internal index identifier
	//            and optional additional information.
	// VcfInfo:   Data specifying a given INFO field
	// VcfFormat: Data specifying a given FORMAT field
	// VcfFilter: Data specifying a given FILTER field
	// VcfStructuredExtra:
	std::vector<std::string>        samples_;
	std::vector<VcfContig>          contigs_;
	// Not written out: used during Import procedure only.
	std::vector<VcfInfo>            info_fields_;
	std::vector<VcfFormat>          format_fields_;
	std::vector<VcfFilter>          filter_fields_;
	std::vector<VcfStructuredExtra> structured_extra_fields_;
	std::vector<VcfExtra>           extra_fields_;

	// Utility members
	//
	// Hash tables allowing the mapping from the unique identifier string
	// (such as contig name) to the relative index offset of that object.
	// This approach requires another layer of indirection when mapping
	// from the index to the actual target. For example:
	//
	// contigs[contigs_map_["chr20"].second] <- maps to the actual target
	//
	// The reverse maps allows the mapping from a unique IDX identifier
	// to the relative index offset of that object. As above, this requires
	// an addition indirect lookup to access the desired object. For example
	// mapping the first occuring FORMAT field to its name:
	//
	// reader->vcf_header_.format_fields_[reader->vcf_header_.format_fields_reverse_map_[container.at(0)->d.fmt[0].id]].id
	//
	// map_type hash tables permits mapping string name -> index offset
	// map_reverse_type hash tables permits mapping integer IDX -> index offset
	map_type samples_map_;
	map_type contigs_map_;
	map_type info_fields_map_;
	map_type format_fields_map_;
	map_type filter_fields_map_;
	map_reverse_type contigs_reverse_map_;       // map IDX -> index offset
	map_reverse_type info_fields_reverse_map_;   // map IDX -> index offset
	map_reverse_type format_fields_reverse_map_; // map IDX -> index offset
	map_reverse_type filter_fields_reverse_map_; // map IDX -> index offset
};

}

}



#endif /* TWK_VCFUTILS_H_ */
