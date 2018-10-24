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
	twk_vimport_settings() : c_level(10), block_size(500), input("-"), output("-"){}

	uint8_t c_level;
	uint32_t block_size;
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
