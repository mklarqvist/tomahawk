#ifndef MATH_OUTPUT_STATISTICS_H_
#define MATH_OUTPUT_STATISTICS_H_

#include <limits>

#include "support/type_definitions.h"
#include "tomahawk/two/output_entry.h"

namespace tomahawk{

struct SummaryStatistics{
private:
	typedef SummaryStatistics self_type;

public:
	SummaryStatistics() :
		total(0),
		total_squared(0),
		n_total(0),
		mean(0),
		standard_deviation(0),
		min(std::numeric_limits<double>::max()),
		max(std::numeric_limits<double>::min())
	{

	}

	bool calculate(void){
		if(this->n_total == 0){
			this->mean = 0;
			this->standard_deviation = 0;
			return false;
		}

		this->mean = this->total / this->n_total;

		if(this->n_total > 1){
			this->standard_deviation = sqrt(this->total_squared/this->n_total - (this->total / this->n_total)*(this->total / this->n_total));
		} else this->standard_deviation = 0;

		return true;
	}

	inline double getSigma(void) const{ return(this->standard_deviation); }
	inline double getSigmaSquared(void) const{ return(this->standard_deviation*this->standard_deviation); }

	template <class T> void operator+=(const T& value){
		this->total         += value;
		this->total_squared += value*value;
		this->n_total += 1;
		if(value < this->min) this->min = value;
		if(value > this->max) this->max = value;
	}

	template <class T> void add(const T& value, const double& weight = 1){
		this->total         += value;
		this->total_squared += value*value;
		this->n_total       += weight;
		if(value < this->min) this->min = value;
		if(value > this->max) this->max = value;
	}

	void reset(void){
		this->total         = 0;
		this->total_squared = 0;
		this->n_total       = 0;
		this->mean          = 0;
		this->standard_deviation = 0;
		this->min = std::numeric_limits<double>::max();
		this->max = std::numeric_limits<double>::min();
	}

	// Accessor functions
	inline double getTotal(void) const{ return(this->total); }
	inline double getTotalSquared(void) const{ return(this->total_squared); }
	inline double getCount(void) const{ return(this->n_total); }
	inline double getMean(void) const{ return(this->mean); }
	inline double getStandardDeviation(void) const{ return(this->standard_deviation); }
	inline double getMin(void) const{ return(this->min); }
	inline double getMax(void) const{ return(this->max); }

	friend std::ostream& operator<<(std::ostream& stream, self_type& self){
		self.calculate();
		stream << self.n_total << "\t" << self.mean << "\t" << self.standard_deviation << "\t" << self.min << "\t" << self.max;
		return(stream);
	}

public:
	double total;
	double total_squared;
	double n_total;
	double mean;
	double standard_deviation;
	double min;
	double max;
};

struct SummaryStatisticsObject{
public:
	typedef SummaryStatisticsObject self_type;
	typedef SummaryStatistics       value_type;

public:
	SummaryStatisticsObject()  = default;
	~SummaryStatisticsObject() = default;

	void operator+=(const io::OutputEntry& entry){
		this->p1 += entry.p1;
		this->p2 += entry.p2;
		this->q1 += entry.q1;
		this->q2 += entry.q2;
		this->R  += entry.R;
		this->R2 += entry.R2;
		this->P  += entry.P;
		this->D  += entry.D;
		this->Dprime += entry.Dprime;
	}

	void calculate(void){
		this->R.calculate();
		this->R2.calculate();
		this->P.calculate();
		this->p1.calculate();
		this->p2.calculate();
		this->q1.calculate();
		this->q2.calculate();
		this->D.calculate();
		this->Dprime.calculate();
	}

	void reset(void){
		this->R.reset();
		this->R2.reset();
		this->P.reset();
		this->p1.reset();
		this->p2.reset();
		this->q1.reset();
		this->q2.reset();
		this->D.reset();
		this->Dprime.reset();
	}

	std::ostream& print(std::ostream& stream, const std::string& family_type, const std::string& contigIDA, const int64_t& position, const std::string& contigIDB){
		stream << family_type << "\tR\t" << contigIDA << "\t" << position << "\t" << contigIDB << "\t";
		stream << this->R << "\n";
		stream << family_type << "\tR2\t" << contigIDA << "\t" << position << "\t" << contigIDB << "\t";
		stream << this->R2 << "\n";
		stream << family_type << "\tP\t" << contigIDA << "\t" << position << "\t" << contigIDB << "\t";
		stream << this->P << "\n";
		stream << family_type << "\tp1\t" << contigIDA << "\t" << position << "\t" << contigIDB << "\t";
		stream << this->p1 << "\n";
		stream << family_type << "\tp2\t" << contigIDA << "\t" << position << "\t" << contigIDB << "\t";
		stream << this->p2 << "\n";
		stream << family_type << "\tq1\t" << contigIDA << "\t" << position << "\t" << contigIDB << "\t";
		stream << this->q1 << "\n";
		stream << family_type << "\tq2\t" << contigIDA << "\t" << position << "\t" << contigIDB << "\t";
		stream << this->q2 << "\n";
		stream << family_type << "\tD\t" << contigIDA << "\t" << position << "\t" << contigIDB << "\t";
		stream << this->D << "\n";
		stream << family_type << "\tDprime\t" << contigIDA << "\t" << position << "\t" << contigIDB << "\t";
		stream << this->Dprime << "\n";
		return(stream);
	}

public:
	value_type R, R2, P, p1, p2, q1, q2, D, Dprime;
};

}


#endif /* MATH_OUTPUT_STATISTICS_H_ */
