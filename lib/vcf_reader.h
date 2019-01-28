#ifndef IO_VCF_READER_H_
#define IO_VCF_READER_H_

#include "utility.h"
#include "header.h"
#include "header_internal.h"

// htslib dependencies.
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"

namespace tomahawk{

class VcfReader {
public:
	typedef VcfReader self_type;

public:
	// Singleton design pattern for retrieving a guaranteed unique pointer
	// to a VcfReader. This choice is to prevent inadvertent writes to the
	// target file as another file-handle is accessing it.
	static std::unique_ptr<self_type> FromFile(const std::string& variants_path, uint32_t n_extra_threads = 0){
		htsFile* fp = hts_open(variants_path.c_str(), "r");
		if(n_extra_threads){
			int ret = hts_set_threads(fp, n_extra_threads);
			if(ret < 0){
				std::cerr << utility::timestamp("ERROR") << "Failed to open multiple handles!" << std::endl;
				return nullptr;
			}
		}

		if (fp == nullptr) {
			std::cerr << utility::timestamp("ERROR")  << "Could not open " << variants_path << std::endl;
			return nullptr;
		}

		bcf_hdr_t* header = bcf_hdr_read(fp);
		if (header == nullptr){
			std::cerr << utility::timestamp("ERROR") << "Couldn't parse header for " << fp->fn << std::endl;
			return nullptr;
		}

		return std::unique_ptr<self_type>(new self_type(variants_path, fp, header));
	}

	bool next(const int unpack_level = BCF_UN_ALL){
		if (bcf_read(this->fp_, this->header_, this->bcf1_) < 0) {
			if (bcf1_->errcode) {
				std::cerr << utility::timestamp("ERROR") << "Failed to parse VCF record: " << bcf1_->errcode << std::endl;
				return false;
			} else {
				return false;
			}
		}

		bcf_unpack(this->bcf1_, unpack_level);
		return true;
	}

	bool next(bcf1_t* bcf_entry, const int unpack_level = BCF_UN_ALL){
		if (bcf_read(this->fp_, this->header_, bcf_entry) < 0) {
			if (bcf_entry->errcode) {
				std::cerr << utility::timestamp("ERROR") << "Failed to parse VCF record: " << bcf1_->errcode << std::endl;
				return false;
			} else {
				//std::cerr << utility::timestamp("ERROR") << "Failed to retrieve a htslib bcf1_t record!" << std::endl;
				return false;
			}
		}

		bcf_unpack(bcf_entry, unpack_level);
		return true;
	}

	/**<
	 * Utility function that writes the VcfHeader literals string into
	 * a target output stream. The literals string does NOT contain
	 * sample information or the column header string ("#CHROM...").
	 * @param stream Dst output stream.
	 */
	inline void PrintLiterals(std::ostream& stream) const{ stream << this->vcf_header_.literals_ << std::endl; }

	/**<
	 * Utility function that writes a valid VCF header output string
	 * to the target stream.
	 * @param stream Dst output stream.
	 */
	void PrintVcfHeader(std::ostream& stream) const{
		this->PrintLiterals(stream);
		stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
		if(this->vcf_header_.samples_.size()){
			stream << "\tFORMAT\t";
			stream << this->vcf_header_.samples_[0];
			for(size_t i = 1; i < this->vcf_header_.samples_.size(); ++i)
				stream << "\t" + this->vcf_header_.samples_[i];
		}
		stream << "\n";
	}

private:
	// Private constructor.
	VcfReader(const std::string& variants_path,
              htsFile* fp,
			  bcf_hdr_t* header) :
    fp_(fp),
    header_(header),
    bcf1_(bcf_init())
{
    if (this->header_->nhrec < 1) {
        std::cerr << utility::timestamp("ERROR") << "Empty header, not a valid VCF." << std::endl;
        return;
    }

    // Store the file-format header string
    if (std::string(this->header_->hrec[0]->key) != "fileformat") {
        std::cerr << utility::timestamp("ERROR") << "Not a valid VCF, fileformat needed: " << variants_path << std::endl;
    } else {
    	this->vcf_header_.fileformat_string_ = std::string(this->header_->hrec[0]->key);
    }

    VcfHeaderInternal* hdr = reinterpret_cast<VcfHeaderInternal*>(&vcf_header_);

    // Fill in the contig info for each contig in the VCF header. Directly
    // accesses the low-level C struct because there are no indirection
    // macros/functions by htslib API.
    // BCF_DT_CTG: offset for contig (CTG) information in BCF dictionary (DT).
    const int n_contigs = this->header_->n[BCF_DT_CTG];
    for (int i = 0; i < n_contigs; ++i) {
        const bcf_idpair_t& idPair = this->header_->id[BCF_DT_CTG][i];
        hdr->AddContigInfo(idPair);
    }

    // Iterate through all hrecs (except the first, which was 'fileformat') to
    // populate the rest of the headers.
    for (int i = 1; i < this->header_->nhrec; i++) {
    	const bcf_hrec_t* hrec0 = this->header_->hrec[i];
		switch (hrec0->type) {
		case BCF_HL_CTG:
			// Contigs are populated above, since they store length in the
			// bcf_idinfo_t* structure.
			break;
		case BCF_HL_FLT:
			hdr->AddFilterInfo(hrec0);
			break;
		case BCF_HL_INFO:
			hdr->AddInfo(hrec0);
			break;
		case BCF_HL_FMT:
			hdr->AddFormatInfo(hrec0);
			break;
		case BCF_HL_STR:
			hdr->AddStructuredExtra(hrec0);
			break;
		case BCF_HL_GEN:
			hdr->AddExtra(hrec0);
			break;
		default:
			std::cerr << utility::timestamp("ERROR") << "Unknown hrec0->type: " << hrec0->type << std::endl;
			break;
		}
    }

    // Populate samples info.
    int n_samples = bcf_hdr_nsamples(this->header_);
    for (int i = 0; i < n_samples; i++)
    	hdr->AddSample(std::string(this->header_->samples[i]));

    this->vcf_header_.BuildReverseMaps();

    // Build literal VCF header string for storage.
    kstring_t htxt = {0,0,0};
	bcf_hdr_format(this->header_, 0, &htxt);
	while (htxt.l && htxt.s[htxt.l-1] == '\0') --htxt.l; // kill trailing zeros
	std::string temp = std::string(htxt.s, htxt.l);
	size_t pos = temp.find("#CHROM"); // search for start of column header line
	temp = temp.substr(0, pos);
	this->vcf_header_.literals_ = temp;
	free(htxt.s);
}

public:
	// Public destructor.
	~VcfReader() {
		bcf_destroy(this->bcf1_);
		bcf_hdr_destroy(this->header_);
		hts_close(this->fp_);
	}

public:
	// Contextual representation of vcf header
	VcfHeader vcf_header_;

	// A pointer to the htslib file used to access the VCF data.
	htsFile * fp_;

	// A htslib header data structure obtained by parsing the header of this VCF.
	bcf_hdr_t * header_;

	// htslib representation of a parsed vcf line.
	bcf1_t* bcf1_;
};

}

#endif /* IO_VCF_READER_H_ */
