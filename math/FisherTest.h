#ifndef FISHERTEST_H_
#define FISHERTEST_H_

#include <cmath>
#include <vector>
#include "../TypeDefinitions.h"

namespace Tomahawk {

#define STIRLING_CONSTANT	0.5 * log(2 * M_PI)

class FisherTest{
public:
	FisherTest(const U32 number);
	~FisherTest();

	double StirlingsApproximation(const double v) const;
	double logN(const U32 value) const;
	double fisherTest(const U32 a, const U32 b, const U32 c, const U32 d) const;
	double fisherTestLess(U32 a, U32 b, U32 c, U32 d) const;
	double fisherTestGreater(U32 a, U32 b, U32 c, U32 d) const;
	double chiSquaredTest(U32& a, U32& b, U32& c, U32& d, U32& total) const;

private:
	void Build(void);

private:
	U32 number;
	std::vector<double> logN_values;
};

__attribute__((always_inline))
inline double FisherTest::StirlingsApproximation(const double v) const{
	return(STIRLING_CONSTANT + (v - 0.5) * log(v) - v + 1/(12*v) + 1/(360 * v * v * v));
}

__attribute__((always_inline))
inline double FisherTest::logN(const U32 value) const{
	if(value < this->number)
		return this->logN_values[value];
	else
		return this->StirlingsApproximation(value + 1);
}

__attribute__((always_inline))
inline double FisherTest::fisherTest(const U32 a, const U32 b, const U32 c, const U32 d) const{
	// Rewrite Fisher's 2x2 test in log form
	// return e^x
	return(exp(logN(a+b) + logN(c+d) + logN(a+c) + logN(b+d) - logN(a) - logN(b) - logN(c) - logN(d) - logN(a + b + c + d)));
}

inline double FisherTest::chiSquaredTest(U32& a, U32& b, U32& c, U32& d, U32& total) const{
	const U32 rowSums[2] = {a+b, c+d};
	const U32 colSums[2] = {a+c, b+d};

	return(pow(a - (double)rowSums[0]*colSums[0]/total,2)/(rowSums[0]*colSums[0]/total) +
		   pow(b - (double)rowSums[0]*colSums[1]/total,2)/(rowSums[0]*colSums[1]/total) +
		   pow(c - (double)rowSums[1]*colSums[0]/total,2)/(rowSums[1]*colSums[0]/total) +
		   pow(d - (double)rowSums[1]*colSums[1]/total,2)/(rowSums[1]*colSums[1]/total));
}

} /* namespace Tomahawk */

#endif /* FISHERTEST_H_ */
