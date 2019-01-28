#ifndef TWK_GENOTYPE_ENCODER_H_
#define TWK_GENOTYPE_ENCODER_H_

#include <cstdint>
#include <cassert>

#include "core.h"

namespace tomahawk {

const static uint8_t TWK_GT_MAP[3]  = {2,0,1};
const static uint8_t TWK_GT_FLIP[3] = {1,0,2}; // usage TWK_GT_FLIP[TWK_GT_MAP[.]]
const static uint8_t TWK_GT_FLIP_NONE[3] = {0,1,2}; // usage TWK_GT_FLIP_NONE[TWK_GT_MAP[.]]
#define TWK_GT_PACK(A,B,MISS) ((TWK_GT_MAP[((A) >> 1)] << ((MISS)+1)) | (TWK_GT_MAP[((B) >> 1)]))
#define TWK_GT_PACK_FLIP(A,B,MISS,FLIP) ((FLIP[TWK_GT_MAP[((A) >> 1)]] << ((MISS)+1)) | (FLIP[TWK_GT_MAP[((B) >> 1)]]))
#define TWK_GT_RLE_PACK(REF,LEN,MISS) (((LEN) << (2+2*(MISS))) | (REF))
#define TWK_GT_LIMIT(T,MISS) ((1L << ((T)-2-2*(MISS))) - 1)

// 0: invariant
// 1: missing threshold
// 2: insufficient samples
// 3: mixed ploidy
// 4: no genotypes
// 5: no fmt
// 6: not biallelic
// 7: not snp
// 8: hardy-weinberg limit
static uint64_t TWK_SITES_FILTERED[9];
const static std::string TWK_SITES_FILTERED_NAMES[9] = {"Invariant","Missing threshold","Insufficient samples","Mixed ploidy","No genotypes","No FORMAT","Not biallelic","Not SNP","Hardy-Weinberg threshold"};

struct GenotypeSummary {
public:
	GenotypeSummary(void) :
		base_ploidy(0),
		phase_if_uniform(0),
		mixed_phasing(false),
		invariant(false),
		n_missing(0),
		n_vector_end(0)
	{
		memset(cnt, 0, sizeof(uint64_t)*3);
		memset(hap_cnt, 0, sizeof(uint64_t)*10);
	}

	~GenotypeSummary() = default;

	/**<
	 * Gathers summary statistics for a vector of genotypes
	 * at a given site. Collects information regarding the
	 * number of missing genotypes and count of sentinel
	 * nodes, checks if the phasing is uniform and whether
	 * all the genotypes are identical.
	 * @param n_samples Total number of samples in the input vector.
	 *                  This is equivalent to the samples in the file.
	 * @param fmt       The target htslib format structure.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	template<class T>
	bool Evaluate(const size_t& n_samples, const bcf_fmt_t& fmt){
		if(fmt.p_len == 0) return true;
		assert(fmt.size/fmt.n == sizeof(T));

		// Set the base ploidy. This corresponds to the LARGEST
		// ploidy observed of ANY individual sample at the given
		// locus. If a genotype has a ploidy < base ploidy then
		// it is trailed with the sentinel node symbol to signal
		// that the remainder of the vector is NULL.
		this->base_ploidy = fmt.n;

		// Find first phase
		this->phase_if_uniform = 0;
		// Iterate over genotypes to find the first valid phase
		// continuing the search if the current value is the
		// sentinel node symbol.
		int j = 0;
		for(uint32_t i = 0; i < n_samples; ++i){
			if(VcfGenotype<T>::IsMissing(fmt.p[j]) == true
			   || VcfType<T>::IsVectorEnd(fmt.p[j]) == true
			   || VcfGenotype<T>::IsMissing(fmt.p[j+1]) == true
			   || VcfType<T>::IsVectorEnd(fmt.p[j+1]) == true)
				j += fmt.n;
			else {
				this->phase_if_uniform = fmt.p[j+1] & 1;
				break;
			}
		}
		//const uint32_t ref_phase_pos = j;
		j = 0;

		// Iterate over genotypes to compute summary statistics
		// regarding missingness, number of special sentinel
		// symbols and assess uniformity of phasing.
		for(uint32_t i = 0; i < n_samples; ++i){
			if(VcfGenotype<T>::IsMissing(fmt.p[j+1]) == false
			   && VcfType<int8_t>::IsVectorEnd(fmt.p[j+1]) == false
			   && (fmt.p[j+1] & 1) != this->phase_if_uniform)
			{
				//std::cerr << "is mixed=" << (int)(fmt.p[j+1] & 1) << "!=" << (int)this->phase_if_uniform << std::endl;
				this->mixed_phasing = true;
			}

			// Iterate over the number of chromosomes / individual
			++hap_cnt[(TWK_GT_MAP[fmt.p[j] >> 1] << 2) | (TWK_GT_MAP[fmt.p[j+1] >> 1])];
			for(int k = 0; k < fmt.n; ++k, ++j){
				assert(j < fmt.p_len);

				++cnt[TWK_GT_MAP[fmt.p[j] >> 1]];
				this->n_missing    += VcfGenotype<T>::IsMissing(fmt.p[j]);
				this->n_vector_end += VcfType<T>::IsVectorEnd(fmt.p[j]);
			}
		}
		assert(j == n_samples*2);

		return true;
	}

	inline bool isBaseDiploid(void) const{ return(this->base_ploidy == 2); }

public:
	uint8_t  base_ploidy;
	bool     phase_if_uniform;
	bool     mixed_phasing;
	bool     invariant;
	uint64_t n_missing;
	uint64_t n_vector_end;
	uint64_t cnt[3]; // ref, alt, miss
	uint64_t hap_cnt[10]; // largest = 2->2 = 1010b = 9
};

class GenotypeEncoder {
public:
	//GenotypeEncoder()
	struct GenotypeHelper {
		uint8_t ptype;
		uint32_t cnt;
	};

	static GenotypeHelper AssessGenotypes(const bcf1_t* rec, const bool missing = false){
		assert(rec->n_fmt != 0);
		assert(rec->d.fmt[0].n == 2);
		const uint8_t* data = rec->d.fmt[0].p;

		uint32_t len[3] = {1,1,1};
		uint32_t cnt[3] = {0,0,0};
		uint32_t cum[3] = {0,0,0};

		uint8_t ref = TWK_GT_PACK(data[0],data[1],missing);
		for(int i = 2; i < rec->n_sample*rec->d.fmt[0].n; i+=2){
			uint8_t cur = TWK_GT_PACK(data[i],data[i+1],missing);

			if(ref != cur){
				//std::cerr << " " << length << "|" << (int)ref;
				ref = cur;
				for(int j = 0; j < 3; ++j){ ++cnt[j]; cum[j] += len[j]; }
				memset(len, 0, sizeof(uint32_t)*3);
			}

			if(len[0] == TWK_GT_LIMIT(8,missing)){
				++cnt[0];
				cum[0] += len[0];
				len[0]  = 0;
			}

			if(len[1] == TWK_GT_LIMIT(16,missing)){
				++cnt[1];
				cum[1] += len[1];
				len[1]  = 0;
			}

			if(len[2] == TWK_GT_LIMIT(32,missing)){
				++cnt[2];
				cum[2] += len[2];
				len[2]  = 0;
			}

			for(int j = 0; j < 3; ++j) ++len[j];
		}
		for(int j = 0; j < 3; ++j){ ++cnt[j]; cum[j] += len[j]; assert(cum[j] == rec->n_sample); }
		//std::cerr << "cnt="  << cnt[0] << "," << cnt[1] << "," << cnt[2] << std::endl;
		//std::cerr << "cost=" << cnt[0]*1 << "," << cnt[1]*2 << "," << cnt[2]*4 << std::endl;

		const uint8_t cost_factor[3] = {1,2,4};
		uint8_t min_type = 0; uint32_t min_cost = cnt[0]; uint32_t min_cnt = cnt[0];
		for(int j = 1; j < 3; ++j){
			if(cnt[j]*cost_factor[j] < min_cost){
				min_cost = cnt[j]*cost_factor[j];
				min_type = j;
				min_cnt  = cnt[j];
			}
		}
		//std::cerr << "winner=" << (int)min_type << ":" << min_cost << std::endl;
		GenotypeHelper ret; ret.ptype = min_type; ret.cnt = min_cnt;

		return(ret);
	}

	static bool Encode(const bcf1_t* rec, twk1_t& twk, const twk_vimport_settings& settings){
		GenotypeSummary gt;
		// Do not support mixed ploidy
		if(gt.n_vector_end) return false;
		gt.Evaluate<int8_t>(rec->n_sample,rec->d.fmt[0]);
		const GenotypeHelper ret = GenotypeEncoder::AssessGenotypes(rec,gt.n_missing != 0);

		// Ascertain that there is sufficient amount of samples to reliably calculate LD.
		//uint64_t total = gt.cnt[0] + gt.cnt[1];
		uint64_t total_hap = gt.hap_cnt[0] + gt.hap_cnt[1] + gt.hap_cnt[4] + gt.hap_cnt[5];
		if(total_hap < settings.threshold_miss*rec->n_sample){
			//std::cerr << "filter miss threshold=" << total_hap << "<" << settings.threshold_miss*rec->n_sample << std::endl;
			++TWK_SITES_FILTERED[1];
			return false;
		}

		if(total_hap < 5){
			//std::cerr << "not enough samples=" << rec->n_sample << " with " << gt.n_missing << " miss -> available=" << 2*rec->n_sample-gt.n_missing << "/" << total << "/" << total_hap << std::endl;
			++TWK_SITES_FILTERED[2];
			return false;
		}

		if(gt.n_vector_end){
			++TWK_SITES_FILTERED[3];
			//std::cerr << "twk do not support mixed ploidy sites" << std::endl;
			return false;
		}

		if(gt.hap_cnt[0] == total_hap || gt.hap_cnt[1] == total_hap ||
		   gt.hap_cnt[4] == total_hap || gt.hap_cnt[5] == total_hap)
		{
			//std::cerr << "site is invariant=" << gt.cnt[0] << "," << gt.cnt[1] << " and " << gt.hap_cnt[0] << "," << gt.hap_cnt[1] << "," << gt.hap_cnt[4] << "," << gt.hap_cnt[5] << std::endl;
			//return false;
			if(settings.remove_univariate){
				++TWK_SITES_FILTERED[0];
				//std::cerr << "removing ivariant site=" << rec->rid << ":" << rec->pos+1 << "..." << std::endl;
				return false;
			}
		}

		bool flip_allele = false;
		if(gt.cnt[1] > gt.cnt[0]){
			//std::cerr << "site is flipped=" << gt.cnt[0] << "," << gt.cnt[1] << std::endl;
			flip_allele = true;
			if(settings.flip_major_minor) twk.gt_flipped = true;
			//if(settings.flip_major_minor) std::cerr << "will flip" << std::endl;
		}

		// If mixed phasing then set to unphased
		if(gt.mixed_phasing) twk.gt_phase = 0;
		else twk.gt_phase = gt.phase_if_uniform;

		if(gt.n_missing){
			//std::cerr << "setting missing for " << twk.pos << std::endl;
			twk.gt_missing = true;
		}

		twk.n_hom = gt.hap_cnt[5];
		twk.n_het = gt.hap_cnt[1] + gt.hap_cnt[4];

		switch(ret.ptype){
		case(0): return GenotypeEncoder::Encode_<uint8_t> (rec, ret.cnt, twk, flip_allele && settings.flip_major_minor, gt.n_missing != 0);
		case(1): return GenotypeEncoder::Encode_<uint16_t>(rec, ret.cnt, twk, flip_allele && settings.flip_major_minor, gt.n_missing != 0);
		case(2): return GenotypeEncoder::Encode_<uint32_t>(rec, ret.cnt, twk, flip_allele && settings.flip_major_minor, gt.n_missing != 0);
		}
		return false;
	}

	/**<
	 * Internal encoding function for genotypes. Run-length encodes genotypes in
	 * a variant-centric fashion while constraining their length to some unified
	 * primitive type (word width).
	 * @param rec
	 * @param cnt
	 * @param twk
	 * @param flip_alleles
	 * @param missing
	 * @return
	 */
	template <class int_t>
	static bool Encode_(const bcf1_t* rec,
	                    const uint32_t cnt,
	                    twk1_t& twk,
						const bool flip_alleles,
	                    const bool missing = false)
	{
		assert(rec->n_fmt != 0);
		assert(rec->d.fmt[0].n == 2);
		const uint8_t* data = rec->d.fmt[0].p;
		const uint32_t limit = TWK_GT_LIMIT(sizeof(int_t)*8,missing);
		twk.gt = new twk1_igt_t<int_t>; // new gt container
		twk1_igt_t<int_t>* gt = reinterpret_cast<twk1_igt_t<int_t>*>(twk.gt);
		gt->data = new int_t[cnt]; // new gt data
		gt->n = cnt; // set number runs
		gt->miss = missing;
		twk.gt_ptype = sizeof(int_t); // set ptype

		const uint8_t* flip_map = (flip_alleles ? TWK_GT_FLIP : TWK_GT_FLIP_NONE);

		// count for an and af
		uint32_t ac[3]; memset(ac,0,sizeof(uint32_t)*3);
		++ac[flip_map[TWK_GT_MAP[(data[0] >> 1)]]];
		++ac[flip_map[TWK_GT_MAP[(data[1] >> 1)]]];

		uint32_t len = 1, cumsum = 0, icnt = 0;
		uint8_t  ref = TWK_GT_PACK_FLIP(data[0],data[1],missing,flip_map);

		for(int i = 2; i < rec->n_sample*rec->d.fmt[0].n; i+=2){
			uint8_t cur = TWK_GT_PACK_FLIP(data[i],data[i+1],missing,flip_map);

			++ac[flip_map[TWK_GT_MAP[(data[i] >> 1)]]];
			++ac[flip_map[TWK_GT_MAP[(data[i+1] >> 1)]]];

			if(ref != cur){
				int_t val = TWK_GT_RLE_PACK(ref,len,missing);
				assert((val >> (2+2*missing)) == len);
				gt->data[icnt++] = val;
				ref = cur;
				cumsum += len;
				len = 0;
			}

			if(len == limit){
				gt->data[icnt++] = TWK_GT_RLE_PACK(ref,len,missing);
				cumsum += len;
				len = 0;
			}
			++len;
		}
		// Add last
		int_t val = TWK_GT_RLE_PACK(ref,len,missing);
		assert((val >> (2+2*missing)) == len);
		gt->data[icnt++] = val;
		cumsum += len;
		assert(cumsum == rec->n_sample);
		assert(icnt == cnt);

		uint32_t n_tot_ac = ac[0] + ac[1] + ac[2];
		assert(n_tot_ac == rec->n_sample*2);
		twk.ac = ac[1];
		twk.an = ac[2];
		//std::cerr << "hets=" << hets << std::endl;
		//std::cerr << twk.pos << ":" << twk.n_het << "," << twk.n_hom << " and " << twk.ac << "," << twk.an << std::endl;


		return true;
	}
};


}


#endif /* GENOTYPE_ENCODER_H_ */
