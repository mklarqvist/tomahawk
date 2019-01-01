#include "utility.h"
#include "header.h"

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

VcfHeader::VcfHeader(const VcfHeader& other) :
	fileformat_string_(other.fileformat_string_),
	literals_(other.literals_),
	samples_(other.samples_),
	contigs_(other.contigs_),
	info_fields_(other.info_fields_),
	format_fields_(other.format_fields_),
	filter_fields_(other.filter_fields_),
	structured_extra_fields_(other.structured_extra_fields_),
	extra_fields_(other.extra_fields_)
{
	this->BuildMaps();
	this->BuildReverseMaps();
}

VcfContig* VcfHeader::GetContig(const std::string& name) {
	map_type::const_iterator it = this->contigs_map_.find(name);
	if(it != this->contigs_map_.end()) return(&this->contigs_[it->second]);
	return(nullptr);
}

VcfContig* VcfHeader::GetContig(const int& idx) {
	map_reverse_type::const_iterator it = this->contigs_reverse_map_.find(idx);
	if(it != this->contigs_reverse_map_.end()) return(&this->contigs_[it->second]);
	return(nullptr);
}

VcfInfo* VcfHeader::GetInfo(const std::string& name) {
	map_type::const_iterator it = this->info_fields_map_.find(name);
	if(it != this->info_fields_map_.end()) return(&this->info_fields_[it->second]);
	return(nullptr);
}

VcfInfo* VcfHeader::GetInfo(const int& idx) {
	map_reverse_type::const_iterator it = this->info_fields_reverse_map_.find(idx);
	if(it != this->info_fields_reverse_map_.end()) return(&this->info_fields_[it->second]);
	return(nullptr);
}

VcfFormat* VcfHeader::GetFormat(const std::string& name) {
	map_type::const_iterator it = this->format_fields_map_.find(name);
	if(it != this->format_fields_map_.end()) return(&this->format_fields_[it->second]);
	return(nullptr);
}

VcfFormat* VcfHeader::GetFormat(const int& idx) {
	map_reverse_type::const_iterator it = this->format_fields_reverse_map_.find(idx);
	if(it != this->format_fields_reverse_map_.end()) return(&this->format_fields_[it->second]);
	return(nullptr);
}

VcfFilter* VcfHeader::GetFilter(const std::string& name) {
	map_type::const_iterator it = this->filter_fields_map_.find(name);
	if(it != this->filter_fields_map_.end()) return(&this->filter_fields_[it->second]);
	return(nullptr);
}

VcfFilter* VcfHeader::GetFilter(const int& idx) {
	map_reverse_type::const_iterator it = this->filter_fields_reverse_map_.find(idx);
	if(it != this->filter_fields_reverse_map_.end()) return(&this->filter_fields_[it->second]);
	return(nullptr);
}

std::string* VcfHeader::GetSample(const std::string& name) {
	map_type::const_iterator it = this->samples_map_.find(name);
	if(it != this->samples_map_.end()) return(&this->samples_[it->second]);
	return(nullptr);
}

const VcfContig* VcfHeader::GetContig(const std::string& name) const {
	map_type::const_iterator it = this->contigs_map_.find(name);
	if(it != this->contigs_map_.end()) return(&this->contigs_[it->second]);
	return(nullptr);
}

const VcfContig* VcfHeader::GetContig(const int& idx) const {
	map_reverse_type::const_iterator it = this->contigs_reverse_map_.find(idx);
	if(it != this->contigs_reverse_map_.end()) return(&this->contigs_[it->second]);
	return(nullptr);
}

const VcfInfo* VcfHeader::GetInfo(const std::string& name) const {
	map_type::const_iterator it = this->info_fields_map_.find(name);
	if(it != this->info_fields_map_.end()) return(&this->info_fields_[it->second]);
	return(nullptr);
}

const VcfInfo* VcfHeader::GetInfo(const int& idx) const {
	map_reverse_type::const_iterator it = this->info_fields_reverse_map_.find(idx);
	if(it != this->info_fields_reverse_map_.end()) return(&this->info_fields_[it->second]);
	return(nullptr);
}

const VcfFormat* VcfHeader::GetFormat(const std::string& name) const {
	map_type::const_iterator it = this->format_fields_map_.find(name);
	if(it != this->format_fields_map_.end()) return(&this->format_fields_[it->second]);
	return(nullptr);
}

const VcfFormat* VcfHeader::GetFormat(const int& idx) const {
	map_reverse_type::const_iterator it = this->format_fields_reverse_map_.find(idx);
	if(it != this->format_fields_reverse_map_.end()) return(&this->format_fields_[it->second]);
	return(nullptr);
}

const VcfFilter* VcfHeader::GetFilter(const std::string& name) const {
	map_type::const_iterator it = this->filter_fields_map_.find(name);
	if(it != this->filter_fields_map_.end()) return(&this->filter_fields_[it->second]);
	return(nullptr);
}

const VcfFilter* VcfHeader::GetFilter(const int& idx) const {
	map_reverse_type::const_iterator it = this->filter_fields_reverse_map_.find(idx);
	if(it != this->filter_fields_reverse_map_.end()) return(&this->filter_fields_[it->second]);
	return(nullptr);
}

const std::string* VcfHeader::GetSample(const std::string& name) const {
	map_type::const_iterator it = this->samples_map_.find(name);
	if(it != this->samples_map_.end()) return(&this->samples_[it->second]);
	return(nullptr);
}

bool VcfHeader::BuildReverseMaps(void){
	this->contigs_reverse_map_.clear();
	this->info_fields_reverse_map_.clear();
	this->format_fields_reverse_map_.clear();
	this->filter_fields_reverse_map_.clear();

	for(uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_reverse_map_[this->contigs_[i].idx] = i;
	for(uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_reverse_map_[this->info_fields_[i].idx] = i;
	for(uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_reverse_map_[this->format_fields_[i].idx] = i;
	for(uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_reverse_map_[this->filter_fields_[i].idx] = i;

	return true;
}

bool VcfHeader::BuildMaps(void){
	this->info_fields_map_.clear();
	this->format_fields_map_.clear();
	this->filter_fields_map_.clear();
	this->contigs_map_.clear();

	for(uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_map_[this->contigs_[i].name] = i;
	for(uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_map_[this->info_fields_[i].id] = i;
	for(uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_map_[this->format_fields_[i].id] = i;
	for(uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_map_[this->filter_fields_[i].id] = i;

	return true;
}

twk_buffer_t& operator<<(twk_buffer_t& buffer, const VcfHeader& self){
	SerializeString(self.fileformat_string_, buffer);
	SerializeString(self.literals_, buffer);

	// Samples
	const uint32_t n_samples = self.samples_.size();
	SerializePrimitive(n_samples, buffer);
	for(int i = 0; i < n_samples; ++i) SerializeString(self.samples_[i], buffer);

	// Contigs
	const uint32_t n_contigs = self.contigs_.size();
	SerializePrimitive(n_contigs, buffer);
	for(int i = 0; i < n_contigs; ++i) buffer << self.contigs_[i];

	return(buffer);
}

twk_buffer_t& operator>>(twk_buffer_t& buffer, VcfHeader& self){
	DeserializeString(self.fileformat_string_, buffer);
	DeserializeString(self.literals_, buffer);

	// Samples
	uint32_t n_samples = 0;
	DeserializePrimitive(n_samples, buffer);
	self.samples_.resize(n_samples);
	for(int i = 0; i < n_samples; ++i) DeserializeString(self.samples_[i], buffer);

	// Contigs
	uint32_t n_contigs = 0;
	DeserializePrimitive(n_contigs, buffer);
	self.contigs_.resize(n_contigs);
	for(int i = 0; i < n_contigs; ++i) buffer >> self.contigs_[i];

	self.BuildMaps();
	self.BuildReverseMaps();

	return(buffer);
}

}
