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
#ifndef TWK_IMPORTER_H_
#define TWK_IMPORTER_H_

#include <cstdint>
#include <string>

namespace tomahawk {

/**<
 * Struct encapsulating the settable parameters for the `twk_variant_importer`
 * class below. If you want to write to stdout then set output to '-'. If reading
 */
struct twk_vimport_settings {
	twk_vimport_settings() : remove_univariate(true), flip_major_minor(false), c_level(1), block_size(500), threshold_miss(0.9), hwe(0), input("-"), output("-"){}

	bool remove_univariate, flip_major_minor;
	uint8_t c_level;
	uint32_t block_size;
	float threshold_miss;
	double hwe;
	std::string input, output;
};

/**<
 * Basic class for importing htslib-compatible files into the tomahawk file format.
 * Importing require the use of the `twk_vimport_settings` struct for providing
 * parameters.
 */
class twk_variant_importer {
public:
	bool Import(twk_vimport_settings& settings);
	bool Import(void);

public:
	twk_vimport_settings settings;
};

}

#endif /* TWK_IMPORTER_H_ */
