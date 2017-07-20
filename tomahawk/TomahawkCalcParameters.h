#ifndef TOMAHAWK_TOMAHAWKCALCPARAMETERS_H_
#define TOMAHAWK_TOMAHAWKCALCPARAMETERS_H_

#include <thread>
#include "../io/BasicWriters.h"

namespace Tomahawk{

#define CALC_DEFAULT_MINR2 0.1
#define CALC_DEFAULT_MAXR2 1
#define CALC_DEFAULT_MINP  1e-4
#define CALC_DEFAULT_MINALLELES 5
#define CALC_DEFAULT_MAXALLELES std::numeric_limits<int64_t>::max()

struct TomahawkCalcParameters{
	typedef TomahawkCalcParameters self_type;
	typedef IO::GenericWriterInterace writer_type;
	enum force_method {none, phasedFunction, unphasedFunction};

	TomahawkCalcParameters() :
		n_threads(std::thread::hardware_concurrency() > 0 ? std::thread::hardware_concurrency() : 1),
		n_chunks(1),
		chunk_selected(0),
		R2_min(CALC_DEFAULT_MINR2),
		R2_max(CALC_DEFAULT_MAXR2),
		P_threshold(CALC_DEFAULT_MINP),
		minimum_alleles(CALC_DEFAULT_MINALLELES),
		maximum_alleles(CALC_DEFAULT_MAXALLELES),
		compression_type(writer_type::compression::binary),
		output_stream_type(writer_type::type::cout),
		force(force_method::none),
		detailed_progress(false)
	{}

	TomahawkCalcParameters(const self_type& other):
		n_threads(other.n_threads),
		n_chunks(other.n_chunks),
		chunk_selected(other.chunk_selected),
		R2_min(other.R2_min),
		R2_max(other.R2_max),
		P_threshold(other.P_threshold),
		minimum_alleles(other.minimum_alleles),
		maximum_alleles(other.maximum_alleles),
		compression_type(other.compression_type),
		output_stream_type(other.output_stream_type),
		force(other.force),
		detailed_progress(other.detailed_progress)
	{}

	~TomahawkCalcParameters(){}

	bool Validate(void);

	friend std::ostream& operator<<(std::ostream& os, const self_type& p){
		// Todo: extend to dump all implicit parameters and store in TWO header output
		os << Helpers::timestamp("CALC", "PARAMETERS") << "R-squared (" << p.R2_min << '-' << p.R2_max << "), P < " << p.P_threshold << ", minor haplotype frequency > " << p.minimum_alleles;
		return(os);
	}

	S32 n_threads;
	S32 n_chunks;
	S32 chunk_selected;
	double R2_min;
	double R2_max;
	double P_threshold;
	int64_t minimum_alleles;
	int64_t maximum_alleles;
	writer_type::compression compression_type;
	writer_type::type output_stream_type;
	force_method force;
	bool detailed_progress;
};

}

#endif /* TOMAHAWK_TOMAHAWKCALCPARAMETERS_H_ */
