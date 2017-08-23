#include "utility.h"

int concat(int argc, char** argv){
	// 1: Input list of TWO files: check files for header and EOF marker
	// Check validity before beginning copy: do not want to concat
	// several large files and then fail hours later because a file was
	// truncated or erroneous
	// Also have to check that sample and contig information is identical
	// or at least mappable to each other
	// 2: Open first and hard copy from 0->size-EOF
	// 3: Next file: open and copy from header->size-eof
	// 4: Last file: open and copy from header->eof

	return(0);
}
