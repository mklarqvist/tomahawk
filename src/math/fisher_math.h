#ifndef FISHERTEST_H_
#define FISHERTEST_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <cerrno>
#include "../support/type_definitions.h"

namespace Tomahawk {
namespace Algorithm {

#define KF_GAMMA_EPS      1e-14
#define KF_TINY           1e-290
#define FISHER_TINY       1e-279
#define STIRLING_CONSTANT 0.5 * log(2 * M_PI)

class FisherMath{

public:
	FisherMath(const U32 number);
	~FisherMath();
	void Build(void);
	double kf_lgamma(double z) const;
	double _kf_gammap(double s, double z) const;
	double _kf_gammaq(double s, double z) const;

	inline double kf_gammap(double s, double z) const{return z <= 1. || z < s? _kf_gammap(s, z) : 1. - _kf_gammaq(s, z);}
	inline double kf_gammaq(double s, double z) const{return z <= 1. || z < s? 1. - _kf_gammap(s, z) : _kf_gammaq(s, z);}


	__attribute__((always_inline))
	inline double StirlingsApproximation(const double v) const{
		return(STIRLING_CONSTANT + (v - 0.5) * log(v) - v + 1/(12*v) + 1/(360 * v * v * v));
	}

	__attribute__((always_inline))
	inline double logN(const S32& value) const{
		return(StirlingsApproximation((double)value+1));
	}

	__attribute__((always_inline))
	inline double fisherTest(const S32& a, const S32& b, const S32& c, const S32 d) const{
		// Rewrite Fisher's 2x2 test in log form
		// return e^x
		return(exp(logN(a+b) + logN(c+d) + logN(a+c) + logN(b+d) - logN(a) - logN(b) - logN(c) - logN(d) - logN(a + b + c + d)));
	}


	double fisherTestLess(S32 a, S32 b, S32 c, S32 d) const;
	double fisherTestGreater(S32 a, S32 b, S32 c, S32 d) const;
	double chisqr(const S32& Dof, const double& Cv) const;

	template <class T>
	double chiSquaredTest(T& a, T& b, T& c, T& d) const;

private:
	U32 number;
	std::vector<double> logN_values;
};




template <class T>
double FisherMath::chiSquaredTest(T& a, T& b, T& c, T& d) const{
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

}
} /* namespace Tomahawk */

#endif /* FISHERTEST_H_ */
