#include <string>

#include "../reader.h"
#include "VCFParser.h"
#include "../../tomahawk/TomahawkImportWriter.h"

namespace Tomahawk {
namespace VCF{

#define DEFAULT_MISSINGNESS_CUTOFF 0.2

VCFParser::VCFParser(std::string inputFile, std::string outputPrefix) :
	outputPrefix(outputPrefix),
	reader_(inputFile)
{}

VCFParser::~VCFParser(){}

bool VCFParser::Build(){
	if(!this->reader_.open()){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Failed to open file..." << std::endl;
		return false;
	}

	if(!this->header_.parse(this->reader_)){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Failed to parse VCF..." << std::endl;
		exit(1);
	}
	if(!this->header_.good()){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Failed to parse VCF (" << this->header_.error_bit << ")..." << std::endl;
		return false;
	}

	if(this->header_.samples == 0){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "No samples detected..." << std::endl;
		return false;
	}

	if(this->header_.samples == 1){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Cannot run " << Tomahawk::Constants::PROGRAM_NAME << " with a single sample..." << std::endl;
		return false;
	}

	// Parse lines
	line_type line(this->header_.size());
	U32* retValue = nullptr;
	U32* prevRetValue = nullptr;
	U32 previous_position = 0;
	U32 previous_position_simple = 0;

	this->reader_.clear();
	this->writer_.setHeader(this->header_);
	// Todo
	if(!this->writer_.Open(this->outputPrefix))
		return false;

	// While there are lines
	while(this->reader_.getLine()){
		// Parse a VCF line
		if(!line.Parse(&this->reader_[0], this->reader_.size())){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not parse..." << std::endl;
			return false;
		}

		// Try to get contig information from header
		if(!this->header_.getContig(std::string(line.CHROM, line.lCHROM), retValue)){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Contig does not exist in header..." << std::endl;
			return false;
		}

		// If there was no previous contig value
		// set it to current value
		if(prevRetValue == nullptr)
			this->writer_.getTotempoleEntry().contigID = *retValue;

		// Switch in chromosome detected
		if((prevRetValue != nullptr) && (prevRetValue != retValue)){
			if(*retValue < *prevRetValue){
				std::cerr << Helpers::timestamp("ERROR", "VCF") << "Contigs are not sorted (" << this->header_[*prevRetValue].name << " > " << this->header_[*retValue].name << ")..." << std::endl;
				exit(1);
			}

			//std::cerr << Helpers::timestamp("DEBUG", "VCF") << "Switch detected: " << this->header_.getContig(*prevRetValue).name << "->" << this->header_.getContig(*retValue).name << "..." << std::endl;
			previous_position = 0;
			previous_position_simple = 0;

			// Get new contig value from header
			// and flush out data
			++this->header_.getContig(*retValue);
			this->writer_.flush();

			// Update index values
			this->writer_.TotempoleSwitch(*retValue, 0);
		}

		if(line.position > this->header_.getContig(*retValue).length){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << this->header_[*retValue].name << ':' << line.position << " > reported max size of contig (" << this->header_[*retValue].length << ")..." << std::endl;
			return false;
		}

		if(line.position < previous_position){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "File is not sorted by coordinates (" << this->header_[*retValue].name << ':' << line.position << " > " << this->header_[*retValue].name << ':' << previous_position << ")..." << std::endl;
			return false;
		}

		const float missing = line.getMissingness(this->header_.samples);
		if(line.position == previous_position && *retValue == *prevRetValue){
			if(!SILENT)
				std::cerr << Helpers::timestamp("WARNING", "VCF") << "Duplicate position (" << this->header_[*retValue].name << ":" << line.position << "): Dropping..." << std::endl;

			goto next;
		}

		// Execute only if the line is simple (biallelic and SNP)
		if(line.IsSimple()){
			if(missing > DEFAULT_MISSINGNESS_CUTOFF){
				if(!SILENT)
					std::cerr << Helpers::timestamp("WARNING", "VCF") << "Large missingness (" << this->header_[*retValue].name << ":" << line.position << ", " << missing*100 << "%).  Dropping..." << std::endl;

				goto next;
			}

			// Flush if output is large
			if(this->writer_.checkSize()){
				++this->header_.getContig(*retValue); // update block count for this contigID
				this->writer_.flush();

				this->writer_.TotempoleSwitch(*retValue, previous_position_simple);
			}
			this->writer_ += line;
		}

		next:
		previous_position = line.position;
		previous_position_simple = line.position;
		prevRetValue = retValue;
		this->reader_.clear();
	} // end while there are vcf lines

	// This only happens if there are no valid entries in the file
	if(retValue == nullptr){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Did not import any variants..." << std::endl;
		return false;
	}

	++this->header_.getContig(*retValue);
	this->writer_.flush();
	//		return false;

	this->writer_.WriteFinal();

	if(this->writer_.GetVariantsWritten() == 0){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Did not import any variants..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Wrote: " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.GetVariantsWritten()))
														 << " variants to " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
														 << " blocks..." << std::endl;

	return true;
}

}
} /* namespace Tomahawk */
