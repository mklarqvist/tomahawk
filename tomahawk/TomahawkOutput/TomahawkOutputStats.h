#ifndef TOMAHAWKOUTPUT_TOMAHAWKOUTPUTSTATS_H_
#define TOMAHAWKOUTPUT_TOMAHAWKOUTPUTSTATS_H_

namespace Tomahawk{
namespace TWO{

struct TomahawkOutputStatsData{
	U64 n_total;
	U64* bins_R2;
};

struct TomahawkOutputStats{
	typedef TomahawkOutputStats self_type;
	typedef TomahawkOutputStatsData data_type;

	data_type within;
	data_type across;
	data_type global;
};

}
}

#endif /* TOMAHAWKOUTPUT_TOMAHAWKOUTPUTSTATS_H_ */
