#include "utility.h"

int stats(int argc, char** argv){
	Tomahawk::IO::TomahawkOutputReader r;
	r.Open("/Users/mk21/Desktop/Projects/Cichlid/Data/run_20170807/cichlid_V1.1__20170807.two");
	r.javelinWeights();

	return(0);
}
