#ifndef TWK_GENOTYPE_ENCODER_H_
#define TWK_GENOTYPE_ENCODER_H_

#include <cstdint>
#include <cassert>

#include "core.h"

namespace tomahawk {

const static uint8_t TWK_GT_MAP[3] = {2,0,1};
#define TWK_GT_PACK(A,B,MISS) ((TWK_GT_MAP[((A) >> 1)] << ((MISS)+1)) | (TWK_GT_MAP[((B) >> 1)]))
#define TWK_GT_RLE_PACK(REF,LEN,MISS) (((LEN) << (2+2*(MISS))) | (REF))
#define TWK_GT_LIMIT(T,MISS) ((1L << ((T)-2-2*(MISS))) - 1)

struct GenotypeSummary {
public:
	GenotypeSummary(void) :
		base_ploidy(0),
		phase_if_uniform(0),
		mixed_phasing(false),
		invariant(false),
		n_missing(0),
		n_vector_end(0)
	{}

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
			   || VcfType<T>::IsVectorEnd(fmt.p[j]) == true)
				j += fmt.n;
			else {
				this->phase_if_uniform = fmt.p[j] & 1;
				break;
			}
		}

		// Iterate over genotypes to compute summary statistics
		// regarding missingness, number of special sentinel
		// symbols and assess uniformity of phasing.
		j = 0;
		for(uint32_t i = 0; i < n_samples; ++i){
			if(VcfGenotype<T>::IsMissing(fmt.p[j]) == false
			   && VcfType<int8_t>::IsVectorEnd(fmt.p[j]) == false
			   && (fmt.p[j] & 1) != this->phase_if_uniform)
			{
				this->mixed_phasing = true;
			}

			// Iterate over the number of chromosomes / individual
			for(int k = 0; k < fmt.n; ++k, ++j){
				assert(j < fmt.p_len);

				this->n_missing    += VcfGenotype<T>::IsMissing(fmt.p[j]);
				this->n_vector_end += VcfType<T>::IsVectorEnd(fmt.p[j]);
			}
		}

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

	static bool Encode(const bcf1_t* rec, twk1_t& twk){
		GenotypeSummary gt;
		// Do not support mixed ploidy
		if(gt.n_vector_end) return false;
		gt.Evaluate<int8_t>(rec->n_sample,rec->d.fmt[0]);
		const GenotypeHelper ret = GenotypeEncoder::AssessGenotypes(rec,gt.n_missing != 0);

		// Ascertain that there is sufficient amount of samples to reliably calculate LD.
		if(rec->n_sample - gt.n_missing < 5){
			std::cerr << "not enough samples=" << rec->n_sample << " with " << gt.n_missing << " miss -> " << rec->n_sample-gt.n_missing << std::endl;
			return false;
		}

		// If mixed phasing then set to unphased
		if(gt.mixed_phasing) twk.gt_phase = 0;
		else twk.gt_phase = gt.phase_if_uniform;

		switch(ret.ptype){
		case(0): return GenotypeEncoder::Encode_<uint8_t> (rec, ret.cnt, twk, gt.n_missing != 0);
		case(1): return GenotypeEncoder::Encode_<uint16_t>(rec, ret.cnt, twk, gt.n_missing != 0);
		case(2): return GenotypeEncoder::Encode_<uint32_t>(rec, ret.cnt, twk, gt.n_missing != 0);
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
	 * @param missing
	 * @return
	 */
	template <class int_t>
	static bool Encode_(const bcf1_t* rec,
	                    const uint32_t cnt,
	                    twk1_t& twk,
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

		// count for an and af
		uint32_t ac[3]; memset(ac,0,sizeof(uint32_t)*3);
		++ac[TWK_GT_MAP[(data[0] >> 1)]];
		++ac[TWK_GT_MAP[(data[1] >> 1)]];

		uint32_t len = 1, cumsum = 0, icnt = 0;
		uint8_t  ref = TWK_GT_PACK(data[0],data[1],missing);
		uint32_t hets = 0;

		hets += ((TWK_GT_MAP[(data[0] >> 1)] == 0 && TWK_GT_MAP[(data[1] >> 1)] == 1) || (TWK_GT_MAP[(data[0] >> 1)] == 1 && TWK_GT_MAP[(data[1] >> 1)] == 0));

		for(int i = 2; i < rec->n_sample*rec->d.fmt[0].n; i+=2){
			uint8_t cur = TWK_GT_PACK(data[i],data[i+1],missing);

			++ac[TWK_GT_MAP[(data[i] >> 1)]];
			++ac[TWK_GT_MAP[(data[i+1] >> 1)]];
			hets += ((TWK_GT_MAP[(data[i] >> 1)] == 0 && TWK_GT_MAP[(data[i+1] >> 1)] == 1) || (TWK_GT_MAP[(data[i] >> 1)] == 1 && TWK_GT_MAP[(data[i+1] >> 1)] == 0));

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
		twk.het = hets;
		//std::cerr << "hets=" << hets << std::endl;

		return true;
	}
};


}


#endif /* GENOTYPE_ENCODER_H_ */
