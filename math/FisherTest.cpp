#include "FisherTest.h"

namespace Tomahawk {
namespace Algorithm{

FisherTest::FisherTest(const U32 number) :
		number(number+1),
		logN_values(number, 0)
{
	this->Build();
}

FisherTest::~FisherTest(){}

void FisherTest::Build(void){
	double factorial = 0;
	this->logN_values[0] = 0;
	for(U32 i = 1; i < this->number; ++i){
		factorial += log(i);
		this->logN_values[i] = factorial;
	}
}

double FisherTest::fisherTestLess(U32 a, U32 b, U32 c, U32 d) const{
	U32 minValue = a;
	if(d < minValue) minValue = d;

	double sum = 0;
	for(U32 i = 0; i <= minValue; ++i){
		sum += this->fisherTest(a, b, c, d);
		--a, ++b, ++c, --d;
	}

	return(sum);
}

double FisherTest::fisherTestGreater(U32 a, U32 b, U32 c, U32 d) const{
	U32 minValue = b;
	if(c < minValue) minValue = c;

	double sum = 0;
	for(U32 i = 0; i <= minValue; ++i){
		sum += this->fisherTest(a, b, c, d);
		++a, --b, --c, ++d;
	}

	return(sum);
}

double FisherTest::chisqr(const S32& Dof, const double& Cv){
	if(Cv < 0 || Dof < 1)
		return 0.0;

	const double K = ((double)Dof) * 0.5;
	const double X = Cv * 0.5;
	if(Dof == 2)
		return exp(-1.0 * X);

   return(this->kf_gammaq(K, X));
}


/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 */
double FisherTest::kf_lgamma(double z){
	double x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

// regularized lower incomplete gamma function, by series expansion
double FisherTest::_kf_gammap(double s, double z){
	double sum, x;
	int k;
	for (k = 1, sum = x = 1.; k < 100; ++k) {
		sum += (x *= z / (s + k));
		if (x / sum < KF_GAMMA_EPS) break;
	}
	return exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
}

// regularized upper incomplete gamma function, by continued fraction
double FisherTest::_kf_gammaq(double s, double z){
	int j;
	double C, D, f;
	f = 1. + z - s; C = f; D = 0.;
	// Modified Lentz's algorithm for computing continued fraction
	// See Numerical Recipes in C, 2nd edition, section 5.2
	for (j = 1; j < 100; ++j) {
		double a = j * (s - j), b = (j<<1) + 1 + z - s, d;
		D = b + a * D;
		if (D < KF_TINY) D = KF_TINY;
		C = b + a / C;
		if (C < KF_TINY) C = KF_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
	}
	return exp(s * log(z) - z - kf_lgamma(s) - log(f));
}

}
} /* namespace Tomahawk */
