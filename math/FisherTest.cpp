#include "FisherTest.h"

namespace Tomahawk {

FisherTest::FisherTest(const U32 number) :
		number(number+1),
		logN_values(number, 0)
{
	this->Build();
}

FisherTest::~FisherTest(){}

void FisherTest::Build(){
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

} /* namespace Tomahawk */
