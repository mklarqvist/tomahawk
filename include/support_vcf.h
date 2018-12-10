/*
Copyright (C) 2016-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TOMAHAWK_SUPPORT_VCF_H_
#define TOMAHAWK_SUPPORT_VCF_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "htslib/vcf.h"
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
static inline int yon_float_is_missing(float f)
{
    union { uint32_t i; float f; } u;
    u.f = f;
    return u.i==TWK_FLOAT_MISSING ? 1 : 0;
}

static inline int yon_float_is_vector_end(float f)
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
  static bool IsMissing(float v)  { return yon_float_is_missing(v); }
  static bool IsVectorEnd(float v){ return yon_float_is_vector_end(v); }
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

struct VcfInfo{
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

struct VcfFormat{
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

struct VcfFilter{
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
struct VcfExtra{
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
struct VcfStructuredExtra{
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

}

#endif /* UTILITY_SUPPORT_VCF_H_ */
