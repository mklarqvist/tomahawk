#ifndef TOMAHAWK_TOMAHAWK_CALC_PARAMETERS_H_
#define TOMAHAWK_TOMAHAWK_CALC_PARAMETERS_H_

#include <thread>

#include "io/basic_writers.h"

namespace tomahawk{

#define CALC_DEFAULT_MINR2      0.1
#define CALC_DEFAULT_MAXR2      1.0
#define CALC_DEFAULT_MINP       1
#define CALC_DEFAULT_MINALLELES 1
#define CALC_DEFAULT_MAXALLELES std::numeric_limits<int64_t>::max()

struct TomahawkCalcParameters{
public:
	typedef TomahawkCalcParameters    self_type;
	typedef io::GenericWriterInterace writer_type;
	enum force_method {none, phasedFunction, unphasedFunction};

public:
	TomahawkCalcParameters() :
		window_mode(false),
		fast_mode(false),
		upper_only(false),
		n_threads(std::thread::hardware_concurrency() > 0 ? std::thread::hardware_concurrency() : 1),
		n_chunks(1),
		chunk_selected(0),
		n_window_bases(0),
		R2_min(CALC_DEFAULT_MINR2),
		R2_max(CALC_DEFAULT_MAXR2),
		P_threshold(CALC_DEFAULT_MINP),
		minimum_sum_alternative_haplotype_count(CALC_DEFAULT_MINALLELES),
		maximum_sum_alternative_haplotype_count(CALC_DEFAULT_MAXALLELES),
		compression_type(writer_type::compression::binary),
		force(force_method::none),
		detailed_progress(false)
	{
	}

	TomahawkCalcParameters(const self_type& other):
		window_mode(other.window_mode),
		fast_mode(other.fast_mode),
		upper_only(other.upper_only),
		n_threads(other.n_threads),
		n_chunks(other.n_chunks),
		chunk_selected(other.chunk_selected),
		n_window_bases(other.n_window_bases),
		R2_min(other.R2_min),
		R2_max(other.R2_max),
		P_threshold(other.P_threshold),
		minimum_sum_alternative_haplotype_count(other.minimum_sum_alternative_haplotype_count),
		maximum_sum_alternative_haplotype_count(other.maximum_sum_alternative_haplotype_count),
		compression_type(other.compression_type),
		force(other.force),
		detailed_progress(other.detailed_progress)
	{
	}

	~TomahawkCalcParameters(){}

	bool Validate(void){
		if(n_threads < 0){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Invalid number of threads..." << std::endl;
			return false;
		}

		if(n_chunks < 0){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Invalid number of partitions..." << std::endl;
			return false;
		}

		if(chunk_selected < 0 || chunk_selected > n_chunks){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Invalid selected partition..." << std::endl;
			return false;
		}

		if(R2_min < 0 || R2_min > 1){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Invalid minimum R-squared cutoff... " << this->R2_min << std::endl;
			return false;
		}

		if(R2_max < 0 || R2_max > 1){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Invalid maximum R-squared cutoff... " << this->R2_max << std::endl;
			return false;
		}

		if(R2_min > R2_max){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Minimum R-squared value > maximum R-squared value..." << std::endl;
			return false;
		}

		if(P_threshold < 0 || P_threshold > 1){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Invalid P-value cutoff..." << std::endl;
			return false;
		}

		if(minimum_sum_alternative_haplotype_count < 0){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Invalid minimum number of alleles..." << std::endl;
			return false;
		}

		if(maximum_sum_alternative_haplotype_count < 0){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Invalid maximum number of alleles..." << std::endl;
			return false;
		}

		if(minimum_sum_alternative_haplotype_count > maximum_sum_alternative_haplotype_count){
			std::cerr << helpers::timestamp("ERROR", "CALC") << "Minimum number of alleles > maximum number of alleles..." << std::endl;
			return false;
		}

		this->R2_min -= constants::ALLOWED_ROUNDING_ERROR;
		this->R2_max += constants::ALLOWED_ROUNDING_ERROR;

		return true;
	}

	std::string getInterpretedString(void) const{
		return(std::string("fast_mode=" + (this->fast_mode ? std::string("TRUE") : std::string("FALSE") ) +
				" window_mode=" + (this->window_mode ? std::string("TRUE window_bases=") + std::to_string(this->n_window_bases) : std::string("FALSE") ) +
				" minR2=" + std::to_string(this->R2_min) + " maxR2=" + std::to_string(this->R2_max) +
				" minP=" + std::to_string(this->P_threshold) +
				" minAHC=" + std::to_string(this->minimum_sum_alternative_haplotype_count) + " maxAHC=" + std::to_string(this->maximum_sum_alternative_haplotype_count) +
				" partStart=" + std::to_string(this->chunk_selected) + " parts="  + std::to_string(this->n_chunks) +
				" threads=" + std::to_string(this->n_threads) + " compression=" + std::to_string(this->compression_type) +
				" force_type=" + std::to_string(this->force)
		));
	}

	friend std::ostream& operator<<(std::ostream& os, const self_type& p){
		// Todo: extend to dump all implicit parameters and store in TWO header output
		os << helpers::timestamp("CALC", "PARAMETERS") << p.getInterpretedString();

		return(os);
	}

public:
	bool    window_mode;
	bool    fast_mode;
	bool    upper_only;
	S32     n_threads;
	S32     n_chunks;
	S32     chunk_selected;
	U32     n_window_bases;
	double  R2_min;
	double  R2_max;
	double  P_threshold;
	int64_t minimum_sum_alternative_haplotype_count;
	int64_t maximum_sum_alternative_haplotype_count;
	writer_type::compression compression_type;
	force_method force;
	bool    detailed_progress;
};

}

#endif /* TOMAHAWK_TOMAHAWK_CALC_PARAMETERS_H_ */
