#ifndef FISHERTEST_H_
#define FISHERTEST_H_

#include <iostream>
#include <cmath>
#include <vector>
#include "../TypeDefinitions.h"

namespace Tomahawk {
namespace Algorithm{

#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290
#define STIRLING_CONSTANT	0.5 * log(2 * M_PI)

class FisherTest{

public:
	FisherTest(const U32 number);
	~FisherTest();

	inline double StirlingsApproximation(const double v) const;
	inline double logN(const U32 value) const;
	inline double fisherTest(const U32 a, const U32 b, const U32 c, const U32 d) const;
	double fisherTestLess(U32 a, U32 b, U32 c, U32 d) const;
	double fisherTestGreater(U32 a, U32 b, U32 c, U32 d) const;
	double chisqr(const S32& Dof, const double& Cv);
	double kf_lgamma(double z);
	double _kf_gammap(double s, double z);
	double _kf_gammaq(double s, double z);
	inline double kf_gammap(double s, double z);
	inline double kf_gammaq(double s, double z);

	template <class T>
	double chiSquaredTest(T& a, T& b, T& c, T& d) const{
		const T rowSums[2] = {a+b, c+d};
		const T colSums[2] = {a+c, b+d};
		const double total = a + b + c + d;
		double adjustValue = 0;
		if(a > 0.5 && b > 0.5 && c > 0.5 && d > 0.5)
			adjustValue = 0.5;

		const double chisq = ((pow(a - (double)rowSums[0]*colSums[0]/total,2)-adjustValue)/(rowSums[0]*colSums[0]/total) +
							  (pow(b - (double)rowSums[0]*colSums[1]/total,2)-adjustValue)/(rowSums[0]*colSums[1]/total) +
							  (pow(c - (double)rowSums[1]*colSums[0]/total,2)-adjustValue)/(rowSums[1]*colSums[0]/total) +
							  (pow(d - (double)rowSums[1]*colSums[1]/total,2)-adjustValue)/(rowSums[1]*colSums[1]/total) );

		if(chisq < 0){
			return ((pow(a - (double)rowSums[0]*colSums[0]/total,2))/(rowSums[0]*colSums[0]/total) +
					(pow(b - (double)rowSums[0]*colSums[1]/total,2))/(rowSums[0]*colSums[1]/total) +
					(pow(c - (double)rowSums[1]*colSums[0]/total,2))/(rowSums[1]*colSums[0]/total) +
					(pow(d - (double)rowSums[1]*colSums[1]/total,2))/(rowSums[1]*colSums[1]/total) );
		}
		return chisq;
	}

private:
	void Build(void);

private:
	U32 number;
	std::vector<double> logN_values;
};

__attribute__((always_inline))
inline double FisherTest::kf_gammap(double s, double z){return z <= 1. || z < s? _kf_gammap(s, z) : 1. - _kf_gammaq(s, z);}

__attribute__((always_inline))
inline double FisherTest::kf_gammaq(double s, double z){return z <= 1. || z < s? 1. - _kf_gammap(s, z) : _kf_gammaq(s, z);}

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

}
} /* namespace Tomahawk */

#endif /* FISHERTEST_H_ */
