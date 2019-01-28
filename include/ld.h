/*
Copyright (C) 2016-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

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
==============================================================================*/
#ifndef TWK_LD_H_
#define TWK_LD_H_

#include <cassert>
#include <thread>

#include "core.h"
#include "twk_reader.h"
#include "writer.h"

namespace tomahawk {

struct twk_ld_balancer; // forward declare balancer

/****************************
*  LD handler
****************************/
class twk_ld {
public:
	twk_ld();
	~twk_ld();

	void operator=(const twk_ld_settings& settings){ this->settings = settings; }

	/**<
	 * Helper function to call Compute subroutine when passing a new settings
	 * object to be used.
	 * @param settings Src settings objects.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool Compute(const twk_ld_settings& settings);
	bool ComputeSingle(const twk_ld_settings& settings, bool verbose = false, bool progress = true);

	/**<
	 * Main subroutine for computing linkage-disequilibrium as contextually
	 * determined given the user-defined parameters.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool Compute();
	bool ComputeSingle(bool verbose = false, bool progress = true);
	bool ComputePerformance();

private:
	class twk_ld_impl;
	twk_ld_settings settings;
	twk_ld_impl* mImpl;
};

}

#endif
