#ifndef TOMAHAWK_TWO_AGGREGATION_PARAMETERS_H_
#define TOMAHAWK_TWO_AGGREGATION_PARAMETERS_H_

namespace tomahawk{
namespace support{

enum TOMAHAWK_AGGREGATION_TARGET{
	TWK_AGGREGATE_COUNT,
	TWK_AGGREGATE_R,
	TWK_AGGREGATE_R2,
	TWK_AGGREGATE_D,
	TWK_AGGREGATE_DPrime,
	TWK_AGGREGATE_P1,
	TWK_AGGREGATE_P2,
	TWK_AGGREGATE_Q1,
	TWK_AGGREGATE_Q2,
	TWK_AGGREGATE_P_VALUE,
	TWK_AGGREGATE_LOG_P_VALUE
};

//count, min, max, mean, sd, sum, sum-squared
enum TOMAHAWK_AGGREGATE_REDUCTION{
	TWK_AGGREGATE_REDUCE_COUNT,
	TWK_AGGREGATE_REDUCE_MEAN,
	TWK_AGGREGATE_REDUCE_MIN,
	TWK_AGGREGATE_REDUCE_MAX,
	TWK_AGGREGATE_REDUCE_SD,
	TWK_AGGREGATE_REDUCE_SUM,
	TWK_AGGREGATE_REDUCE_SUM_SQUARED
};

struct aggregation_parameters{
public:
	aggregation_parameters() :
		scene_x_pixels(0),
		scene_y_pixels(0),
		aggregation_target(TWK_AGGREGATE_R2),
		reduction_target(TWK_AGGREGATE_REDUCE_MEAN)
	{

	}

	~aggregation_parameters() = default;

public:
	U32 scene_x_pixels;
	U32 scene_y_pixels;
	TOMAHAWK_AGGREGATION_TARGET aggregation_target;
	TOMAHAWK_AGGREGATE_REDUCTION reduction_target;
};

}
}



#endif /* TOMAHAWK_TWO_AGGREGATION_PARAMETERS_H_ */
