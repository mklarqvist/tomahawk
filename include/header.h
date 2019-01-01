#ifndef TWK_VCFUTILS_H_
#define TWK_VCFUTILS_H_

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "buffer.h"

#define TWK_BYTE_MISSING   INT8_MIN
#define TWK_BYTE_EOV       (INT8_MIN+1)
#define TWK_SHORT_MISSING  INT16_MIN
#define TWK_SHORT_EOV      (INT16_MIN+1)
#define TWK_INT_MISSING    INT32_MIN
#define TWK_INT_EOV        (INT32_MIN+1)
#define TWK_FLOAT_NAN      0x7FC00000
#define TWK_FLOAT_MISSING  0x7F800001
#define TWK_FLOAT_EOV      0x7F800002
#define TWK_BCF_GT_MISSING 0

// taken from htslib vcf.c and renamed for convenience
static inline int twk_float_is_missing(float f)
{
    union { uint32_t i; float f; } u;
    u.f = f;
    return u.i==TWK_FLOAT_MISSING ? 1 : 0;
}

static inline int twk_float_is_vector_end(float f)
{
    union { uint32_t i; float f; } u;
    u.f = f;
    return u.i==TWK_FLOAT_EOV ? 1 : 0;
}

namespace tomahawk {

//
// -----------------------------------------------------------------------------
// VCF type encoding utilities
template<class T>
struct VcfType {
  // Predicates for checking missing and sentinel entries.  Use these, not ==.
  // Is argument the "missing" value?
  static bool IsMissing(T);
  // Is argument the vector end sentinel value?
  static bool IsVectorEnd(T);
};

// See interface description comment above.
template<>
struct VcfType<int8_t> {
  static bool IsMissing(int8_t v)  { return (v == TWK_BYTE_MISSING); }
  static bool IsVectorEnd(int8_t v){ return (v == TWK_BYTE_EOV); }
};

// See interface description comment above.
template<>
struct VcfType<int16_t> {
  static bool IsMissing(int16_t v)  { return (v == TWK_SHORT_MISSING); }
  static bool IsVectorEnd(int16_t v){ return (v == TWK_SHORT_EOV); }
};

// See interface description comment above.
template<>
struct VcfType<int> {
  static bool IsMissing(int v)  { return (v == TWK_INT_MISSING); }
  static bool IsVectorEnd(int v){ return (v == TWK_INT_EOV); }
};

template <class T>
struct VcfGenotype {
	// Predicates for checking missing and sentinel entries.
	static bool IsMissing(const T& value){ return(value == TWK_BCF_GT_MISSING); }
};

template<>
struct VcfType<float> {
  static bool IsMissing(float v)  { return twk_float_is_missing(v); }
  static bool IsVectorEnd(float v){ return twk_float_is_vector_end(v); }
};


struct VcfContig {
public:
	VcfContig();
	~VcfContig() = default;

	std::string ToVcfString(const bool is_bcf = false) const;

	friend twk_buffer_t& operator<<(twk_buffer_t& buffer, const VcfContig& self){
		SerializePrimitive(self.idx, buffer);
		SerializeString(self.name, buffer);
		SerializeString(self.description, buffer);
		SerializePrimitive(self.n_bases, buffer);
		const uint32_t n_extra = self.extra.size();
		SerializePrimitive(n_extra, buffer);
		for(int i = 0; i < n_extra; ++i){
			SerializeString(self.extra[i].first, buffer);
			SerializeString(self.extra[i].second, buffer);
		}

		return(buffer);
	}

	friend twk_buffer_t& operator>>(twk_buffer_t& buffer, VcfContig& self){
		DeserializePrimitive(self.idx, buffer);
		DeserializeString(self.name, buffer);
		DeserializeString(self.description, buffer);
		DeserializePrimitive(self.n_bases, buffer);
		uint32_t n_extra = 0;
		DeserializePrimitive(n_extra, buffer);
		self.extra.resize(n_extra);
		for(int i = 0; i < n_extra; ++i){
			DeserializeString(self.extra[i].first, buffer);
			DeserializeString(self.extra[i].second, buffer);
		}

		return(buffer);
	}

public:
	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The name of the contig. Canonically this is the first
	// non-whitespace-containing string after the > marker in a FASTA file.
	// For example, the line:
	//      >chr1 more info here
	// has a name of "chr1" and a description of "more info here"
	std::string name;

	// Ideally this record is filled in as described above, but not all FASTA
	// readers capture the description information after the name. Since a
	// description is not required by the FASTA spec, we cannot distinguish cases
	// where a description was not present and where a parser ignored it.
	std::string description;

	// The length of this contig in basepairs.
	int64_t n_bases;

	// Additional information used when reading and writing VCF headers. An
	// example map of key-value extra fields would transform an input line
	// containing 'assembly=B36,taxonomy=x,species="Homo sapiens"' to a map with
	// "assembly" -> "B36", "taxonomy" -> "x", "species" -> "Homo sapiens". We
	// never use this information internally, other than reading it in so we can
	// write the contig out again.
	std::vector< std::pair<std::string, std::string> > extra;
};

struct VcfInfo {
public:
	VcfInfo();
	~VcfInfo() = default;

	std::string ToVcfString(const bool is_bcf = false) const;
	std::string ToVcfString(const uint32_t idx) const;

public:
	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The unique ID of the INFO field. Examples include "MQ0" or "END".
	std::string id;

	// Required. The number of values included with the info field. This should be
	// the string representation of the number, e.g. "1" for a single entry, "2"
	// for a pair of entries, etc. Special cases arise when the number of entries
	// depend on attributes of the Variant or are unknown in advance, and include:
	// "A": The field has one value per alternate allele.
	// "R": The field has one value per allele (including the reference).
	// "G": The field has one value for each possible genotype.
	// ".": The number of values varies, is unknown, or is unbounded.
	std::string number;

	// Required. The type of the INFO field. Valid values are "Integer", "Float",
	// "Flag", "Character", and "String".
	std::string type;

	// Required by VCF. The description of the field.
	std::string description;

	// Optional. The annotation source used to generate the field.
	std::string source;

	// Optional. The version of the annotation source used to generate the field.
	std::string version;
};

struct VcfFormat {
public:
	VcfFormat();
	~VcfFormat() = default;

	std::string ToVcfString(const bool is_bcf = false) const;
	std::string ToVcfString(const uint32_t idx) const;

public:
	// Required. The unique ID of the FORMAT field. Examples include "GT", "PL".
	std::string id;

	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The number of entries expected. See description above in the
	// VcfInfo message.
	std::string number;

	// Required. The type of the field. Valid values are "Integer", "Float",
	// "Character", and "String" (same as INFO except "Flag" is not supported).
	std::string type;

	// Required by VCF. The description of the field.
	std::string description;
};

struct VcfFilter {
public:
	VcfFilter();
	~VcfFilter() = default;

	std::string ToVcfString(const bool is_bcf = false) const;
	std::string ToVcfString(const uint32_t idx) const;

	friend std::ostream& operator<<(std::ostream& stream, const VcfFilter& flt);
	friend std::istream& operator>>(std::istream& stream, VcfFilter& flt);

public:
	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The unique ID of the filter. Examples include "PASS", "RefCall".
	std::string id;

	// Required by VCF. The description of the filter.
	std::string description;
};

// This record type is a catch-all for other types of headers. For example,
// ##pedigreeDB=http://url_of_pedigrees
// The VcfExtra message would represent this with key="pedigreeDB",
// value="http://url_of_pedigrees".
struct VcfExtra {
public:
	VcfExtra() = default;
	VcfExtra(const std::string& key, const std::string& value);
	~VcfExtra() = default;

	std::string ToVcfString(void) const;

	friend std::ostream& operator<<(std::ostream& stream, const VcfExtra& extra);
	friend std::istream& operator>>(std::istream& stream, VcfExtra& extra);

public:
  // Required by VCF. The key of the extra header field. Note that this key does
  // not have to be unique within a VcfHeader.
  std::string key;

  // Required by VCF. The value of the extra header field.
  std::string value;
};

// This record type is a catch-all for other headers containing multiple
// key-value pairs. For example, headers may have META lines that provide
// metadata about the VCF as a whole, e.g.
// ##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
// The VcfStructuredExtra message would represent this with key="META",
// and fields mapping "ID" -> "Assay", "Type" -> "String", etc.
struct VcfStructuredExtra {
public:
	VcfStructuredExtra() = default;
	~VcfStructuredExtra() = default;

	std::string ToVcfString(void) const;

	friend std::ostream& operator<<(std::ostream& stream, const VcfStructuredExtra& extra);
	friend std::istream& operator>>(std::istream& stream, VcfStructuredExtra& extra);

public:
	// Required by VCF. The key of the extra header field. Note that this key does
	// not have to be unique within a VcfHeader.
	std::string key;

	// Required by VCF. The key=value pairs contained in the structure.
	std::vector<VcfExtra> fields;
};

class VcfHeader {
public:
	typedef VcfHeader self_type;
	typedef VcfContig contig_type;
	//typedef bcf_hdr_t hts_vcf_header;
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



#endif /* TWK_VCFUTILS_H_ */
