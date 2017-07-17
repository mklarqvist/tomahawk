#ifndef FISHERTEST_H_
#define FISHERTEST_H_

#include <cmath>
#include <vector>
#include "../TypeDefinitions.h"

namespace Tomahawk {

class FisherTest{
	const double constant1 = 0.5 * log(2 * M_PI);

public:
	FisherTest(const U32 number) : number(number+1), logN_values(number, 0){
		this->Build();
	}
	~FisherTest(){}

	inline double StirlingsApproximation(const double v) const{
		return(constant1 + (v - 0.5) * log(v) - v + 1/(12*v) + 1/(360 * v * v * v));
	}

	inline double logN(const unsigned int value) const{
		if(value < this->number)
			return this->logN_values[value];
		else
			return this->StirlingsApproximation(value + 1);
	}

	double fisherTest(const unsigned int a, const unsigned int b, const unsigned int c, const unsigned int d) const{
	    // Rewrite Fisher's 2x2 test in log form
	    // return e^x
		return(exp(logN(a+b) + logN(c+d) + logN(a+c) + logN(b+d) - logN(a) - logN(b) - logN(c) - logN(d) - logN(a + b + c + d)));
	}

	double fisherTestLess(unsigned int a, unsigned int b, unsigned int c, unsigned int d) const{
	    unsigned int minValue = a;
		if(d < minValue) minValue = d;

		double sum = 0;
		for(unsigned int i = 0; i <= minValue; ++i){
			sum += this->fisherTest(a, b, c, d);
			--a, ++b, ++c, --d;
		}

		return(sum);
	}

	double fisherTestGreater(unsigned int a, unsigned int b, unsigned int c, unsigned int d) const{
	    unsigned int minValue = b;
	    if(c < minValue) minValue = c;

		double sum = 0;
		for(unsigned int i = 0; i <= minValue; ++i){
			sum += this->fisherTest(a, b, c, d);
			++a, --b, --c, ++d;
		}

	    return(sum);
	}

	double chiSquaredTest(U32& a, U32& b, U32& c, U32& d, U32& total) const{
		const U32 rowSums[2] = {a+b, c+d};
		const U32 colSums[2] = {a+c, b+d};

		return(pow(a - (double)rowSums[0]*colSums[0]/total,2)/(rowSums[0]*colSums[0]/total) +
							 pow(b - (double)rowSums[0]*colSums[1]/total,2)/(rowSums[0]*colSums[1]/total) +
							 pow(c - (double)rowSums[1]*colSums[0]/total,2)/(rowSums[1]*colSums[0]/total) +
							 pow(d - (double)rowSums[1]*colSums[1]/total,2)/(rowSums[1]*colSums[1]/total));


	}

private:
	void Build(){
		double factorial = 0;
		this->logN_values[0] = 0;
		for(U32 i = 1; i < this->number; ++i){
			factorial += log(i);
			this->logN_values[i] = factorial;
		}
	}

private:
	U32 number;
	std::vector<double> logN_values;
};

} /* namespace Tomahawk */

#endif /* FISHERTEST_H_ */
