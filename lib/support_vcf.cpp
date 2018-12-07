#include "utility.h"
#include "support_vcf.h"

namespace tomahawk {

VcfContig::VcfContig() : idx(0), n_bases(0){}

std::string VcfContig::ToVcfString(const bool is_bcf) const{
	// Template:
	// ##contig=<ID=GL000241.1,assembly=b37,length=42152>
	std::string ret = "##contig=<ID=" + this->name;
	if(extra.size()){
		ret += "," + this->extra[0].first + "=" + this->extra[0].second;
		for(uint32_t i = 1; i < this->extra.size(); ++i){
			ret += "," + this->extra[i].first + "=" + this->extra[i].second;
		}
	}
	if(this->description.size()) ret += ",Description=" + this->description;
	ret += ",length=" + std::to_string(this->n_bases);
	if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
	ret += ">";
	return(ret);
}

VcfInfo::VcfInfo() : idx(0){}

std::string VcfInfo::ToVcfString(const bool is_bcf) const{
	// Template:
	// ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
	std::string ret = "##INFO=<ID=" + this->id;
	ret += ",Number=" + this->number;
	ret += ",Type=" + this->type;
	ret += ",Description=" + this->description;
	if(this->source.size()) ret += ",Source=" + this->source;
	if(this->source.size()) ret += ",Version=" + this->version;
	if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
	ret += ">";
	return(ret);
}

std::string VcfInfo::ToVcfString(const uint32_t idx) const{
	// Template:
	// ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
	std::string ret = "##INFO=<ID=" + this->id;
	ret += ",Number=" + this->number;
	ret += ",Type=" + this->type;
	ret += ",Description=" + this->description;
	if(this->source.size()) ret += ",Source=" + this->source;
	if(this->source.size()) ret += ",Version=" + this->version;
	ret += ",IDX=" + std::to_string(idx);
	ret += ">";
	return(ret);
}

VcfFormat::VcfFormat() : idx(0){}

std::string VcfFormat::ToVcfString(const bool is_bcf) const{
	// Template:
	// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	std::string ret = "##FORMAT=<ID=" + this->id;
	ret += ",Number=" + this->number;
	ret += ",Type=" + this->type;
	ret += ",Description=" + this->description;
	if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
	ret += ">";
	return(ret);
}

std::string VcfFormat::ToVcfString(const uint32_t idx) const{
	// Template:
	// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	std::string ret = "##FORMAT=<ID=" + this->id;
	ret += ",Number=" + this->number;
	ret += ",Type=" + this->type;
	ret += ",Description=" + this->description;
	ret += ",IDX=" + std::to_string(idx);
	ret += ">";
	return(ret);
}


VcfFilter::VcfFilter() : idx(0){}

std::string VcfFilter::ToVcfString(const bool is_bcf) const{
	// Template:
	// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	std::string ret = "##FILTER=<ID=" + this->id;
	ret += ",Description=" + this->description;
	if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
	ret += ">";
	return(ret);
}

std::string VcfFilter::ToVcfString(const uint32_t idx) const{
	// Template:
	// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	std::string ret = "##FILTER=<ID=" + this->id;
	ret += ",Description=" + this->description;
	ret += ",IDX=" + std::to_string(idx);
	ret += ">";
	return(ret);
}

std::ostream& operator<<(std::ostream& stream, const VcfFilter& flt){
	stream.write((const char*)&flt.idx, sizeof(uint32_t));

	SerializeString(flt.id, stream);
	SerializeString(flt.description, stream);

	return(stream);
}

std::istream& operator>>(std::istream& stream, VcfFilter& flt){
	stream.read((char*)&flt.idx, sizeof(uint32_t));

	DeserializeString(flt.id, stream);
	DeserializeString(flt.description, stream);

	return(stream);
}

VcfExtra::VcfExtra(const std::string& key, const std::string& value) :
	key(key),
	value(value)
{}

std::string VcfExtra::ToVcfString(void) const{
	// Template:
	// ##source=CombineGVCFs
	std::string ret = "##" + this->key + "=" + this->value;
	return(ret);
}

std::ostream& operator<<(std::ostream& stream, const VcfExtra& extra){
	SerializeString(extra.key, stream);
	SerializeString(extra.value, stream);
	return(stream);
}

std::istream& operator>>(std::istream& stream, VcfExtra& extra){
	DeserializeString(extra.key, stream);
	DeserializeString(extra.value, stream);
	return(stream);
}


std::string VcfStructuredExtra::ToVcfString(void) const{
	// Template:
	// ##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
	std::string ret = "##" + this->key + "=<";
	ret += this->fields[0].key + "=" + this->fields[0].value;
	for(uint32_t i = 1; i < this->fields.size(); ++i)
		ret += "," + this->fields[i].key + "=" + this->fields[i].value;
	ret += ">";
	return(ret);
}

std::ostream& operator<<(std::ostream& stream, const VcfStructuredExtra& extra){
	SerializeString(extra.key, stream);
	size_t l_extra = extra.fields.size();
	stream.write((const char*)&l_extra, sizeof(size_t));
	for(uint32_t i = 0; i < extra.fields.size(); ++i)
		stream << extra.fields[i];

	return(stream);
}

std::istream& operator>>(std::istream& stream, VcfStructuredExtra& extra){
	DeserializeString(extra.key, stream);
	size_t l_extra;
	stream.read((char*)&l_extra, sizeof(size_t));
	extra.fields.resize(l_extra);
	for(uint32_t i = 0; i < extra.fields.size(); ++i)
		stream >> extra.fields[i];

	return(stream);
}

}
