#ifndef IMPORTER_H_
#define IMPORTER_H_

#include <cstdint>
#include <string>

namespace tomahawk {

/**<
 * Struct encapsulating the settable parameters for the `twk_variant_importer`
 * class below. If you want to write to stdout then set output to '-'. If reading
 */
struct twk_vimport_settings {
	twk_vimport_settings() : remove_univariate(true), flip_major_minor(false), c_level(1), block_size(500), threshold_miss(0.9), input("-"), output("-"){}

	bool remove_univariate, flip_major_minor;
	uint8_t c_level;
	uint32_t block_size;
	float threshold_miss;
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

#endif /* IMPORTER_H_ */
