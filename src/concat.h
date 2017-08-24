/*
Copyright (C) 2016-2017 Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk21@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/
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
