#include "header_internal.h"

namespace tomahawk {

const bcf_hrec_t* GetPopulatedHrec(const bcf_idpair_t& idPair) {
	for (int i = 0; i < 3; i++) {
		const bcf_hrec_t* hrec = idPair.val->hrec[i];
		if (hrec != nullptr) {
			return hrec;
		}
	}
	std::cerr << "No populated hrec in idPair. Error in htslib." << std::endl;
	return nullptr;
}

// Adds Contig information from the idPair to the ContigInfo object.
void VcfHeaderInternal::AddContigInfo(const bcf_idpair_t& idPair) {
  // ID and length are special-cased in the idPair.
	//std::cerr << "Contig: " << pos_in_fasta << "\t" << std::string(idPair.key) << ": " << idPair.val->info[0] << std::endl;
	VcfContig c;
	c.name = idPair.key;
	c.n_bases = idPair.val->info[0];

  const bcf_hrec_t* hrec0 = GetPopulatedHrec(idPair);
  if (hrec0 != nullptr) {
	for (int j = 0; j < hrec0->nkeys; j++) {
		const std::string current_key(hrec0->keys[j]);
		// Add any non-ID and non-length info to the structured map of additional
		// information.
		if (current_key == "ID" ||
			current_key == "length")
		{
			//continue;
		} else if(current_key == "IDX"){
			c.idx = atoi(hrec0->vals[j]);
		} else {
			c.extra.push_back(std::pair<std::string,std::string>(current_key, std::string(hrec0->vals[j])));
		}
	}
  } else {
	  std::cerr << utility::timestamp("ERROR") << "hrec error" << std::endl;
	  return;
  }

  // Add current contig to map
  if(this->contigs_.size() == 0){
	  this->contigs_.push_back(c);
	  this->contigs_map_[c.name] = 0;
	  return;
  }

  if(this->contigs_map_.find(c.name) == this->contigs_map_.end()){
	  this->contigs_map_[c.name] = this->contigs_.size();
	  this->contigs_.push_back(c);
  } else {
	  std::cerr << utility::timestamp("ERROR") << "Illegal: duplicated contig name" << std::endl;
	  exit(1);
  }
}

// Adds FILTER information from the bcf_hrec_t to the VcfFilterInfo object.
void VcfHeaderInternal::AddFilterInfo(const bcf_hrec_t* hrec) {
  if (hrec->nkeys >= 2 && std::string(hrec->keys[0]) == "ID" &&
	  std::string(hrec->keys[1]) == "Description")
  {
	VcfFilter f;
	f.id = std::string(hrec->vals[0]);
	f.description = std::string(hrec->vals[1]);
	for(int i = 2; i < hrec->nkeys; ++i){
		if(std::string(hrec->keys[i]) == "IDX"){
			f.idx = atoi(hrec->vals[i]);
		}
	}

	// Add current filter field to map.
	if(this->filter_fields_.size() == 0){
		this->filter_fields_.push_back(f);
		this->filter_fields_map_[f.id] = 0;
		return;
	}

	if(this->filter_fields_map_.find(f.id) == this->filter_fields_map_.end()){
		this->filter_fields_map_[f.id] = this->filter_fields_.size();
		this->filter_fields_.push_back(f);
	} else {
		std::cerr << utility::timestamp("ERROR")  << "Illegal: duplicated filter name: " << f.id << std::endl;
		exit(1);
	}

  } else {
	std::cerr << utility::timestamp("ERROR") << "Malformed FILTER field detected in header, leaving this "
				 "filter empty" << std::endl;
  }
}

// Adds INFO information from the bcf_hrec_t to the VcfInfo object.
void VcfHeaderInternal::AddInfo(const bcf_hrec_t* hrec) {
  if (hrec->nkeys >= 4 && std::string(hrec->keys[0]) == "ID" &&
	  std::string(hrec->keys[1]) == "Number" && std::string(hrec->keys[2]) == "Type" &&
	  std::string(hrec->keys[3]) == "Description")
  {
	VcfInfo f;
	f.id = std::string(hrec->vals[0]);
	f.number = std::string(hrec->vals[1]);
	f.type = std::string(hrec->vals[2]);
	f.description = std::string(hrec->vals[3]);
	for (int i = 4; i < hrec->nkeys; i++) {
	  const std::string current_key = std::string(hrec->keys[i]);
	  if (current_key == "Source") {
		f.source = std::string(hrec->vals[i]);
	  } else if (current_key == "Version") {
		f.version = std::string(hrec->vals[i]);
	  } else if (current_key == "IDX") {
		  f.idx = atoi(hrec->vals[i]);
	  }
	}

	// Add current info field to map.
	if(this->info_fields_.size() == 0){
		this->info_fields_.push_back(f);
		this->info_fields_map_[f.id] = 0;
		return;
	}

	if(this->info_fields_map_.find(f.id) == this->info_fields_map_.end()){
		this->info_fields_map_[f.id] = this->info_fields_.size();
		this->info_fields_.push_back(f);
	} else {
		std::cerr << utility::timestamp("ERROR")  << "Illegal: duplicated info name: " << f.id << std::endl;
		exit(1);
	}

  } else {
	std::cerr << utility::timestamp("ERROR") << "Malformed INFO field detected in header, leaving this "
				 "info empty" << std::endl;
  }
}

// Adds FORMAT information from the bcf_hrec_t to the VcfFormatInfo object.
void VcfHeaderInternal::AddFormatInfo(const bcf_hrec_t* hrec) {
  if (hrec->nkeys >= 4 && std::string(hrec->keys[0]) == "ID" &&
	  std::string(hrec->keys[1]) == "Number" && std::string(hrec->keys[2]) == "Type" &&
	  std::string(hrec->keys[3]) == "Description")
  {
	VcfFormat f;
	f.id          = std::string(hrec->vals[0]);
	f.number      = std::string(hrec->vals[1]);
	f.type        = std::string(hrec->vals[2]);
	f.description = std::string(hrec->vals[3]);
	for (int i = 4; i < hrec->nkeys; i++) {
		const std::string current_key = std::string(hrec->keys[i]);
		if (current_key == "IDX") {
			  f.idx = atoi(hrec->vals[i]);
		  }
	}

	// Add current format field to map.
	if(this->format_fields_.size() == 0){
		this->format_fields_.push_back(f);
		this->format_fields_map_[f.id] = 0;
		return;
	}

	if(this->format_fields_map_.find(f.id) == this->format_fields_map_.end()){
		this->format_fields_map_[f.id] = this->format_fields_.size();
		this->format_fields_.push_back(f);
	} else {
		std::cerr << utility::timestamp("ERROR") << "Illegal: duplicated format name: " << f.id << std::endl;
		exit(1);
	}

  } else {
	std::cerr << utility::timestamp("ERROR")  << "Malformed FORMAT field detected in header, leaving this "
					"format empty" << std::endl;
  }
}

// Adds structured information from the bcf_hrec_t to the VcfStructuredExtra.
void VcfHeaderInternal::AddStructuredExtra(const bcf_hrec_t* hrec) {
  VcfStructuredExtra f;
  f.key = std::string(hrec->key);
  for (int i = 0; i < hrec->nkeys; i++)
	  f.fields.push_back(VcfExtra(std::string(hrec->keys[i]), std::string(hrec->vals[i])));

  this->structured_extra_fields_.push_back(f);
}

// Adds unstructured information from the bcf_hrec_t to the VcfExtra object.
void VcfHeaderInternal::AddExtra(const bcf_hrec_t* hrec) {
  VcfExtra f;
  f.key   = std::string(hrec->key);
  f.value = std::string(hrec->value);
  this->extra_fields_.push_back(f);
}

void VcfHeaderInternal::AddSample(const std::string& sample_name) {
	if(this->samples_.size() == 0){
		this->samples_.push_back(sample_name);
		this->samples_map_[sample_name] = 0;
		return;
	}

	if(this->samples_map_.find(sample_name) == this->samples_map_.end()){
		this->samples_map_[sample_name] = this->samples_.size();
		this->samples_.push_back(sample_name);
	} else {
		std::cerr << utility::timestamp("ERROR") << "Illegal: duplicated sample name: " << sample_name << std::endl;
		exit(1);
	}
}

}
