#include <string>

#include "../reader.h"
#include "VCFHeaderConstants.h"
#include "VCFLines.h"
#include "VCFParser.h"
#include "../../tomahawk/TomahawkImportWriter.h"

namespace Tomahawk {
namespace VCF{

#define DEFAULT_MISSINGNESS_CUTOFF 0.2

VCFParser::VCFParser(readerType reader, const std::string outputPrefix) : outputPrefix(outputPrefix), reader_(reader){}
VCFParser::~VCFParser(){}

bool VCFParser::Build(){
	if(!this->GetHeaderLines())
		return false;

	if(!this->header_.BuildContigTable())
		return false;

	if(!this->SampleLine())
		return false;

	if(this->header_.samples_ == 1){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Cannot run " << Tomahawk::Constants::PROGRAM_NAME << " with a single sample..." << std::endl;
		return false;
	}

	// Parse lines
	Tomahawk::VCF::VCFLine line(this->header_.size());
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

		const float missing = line.getMissingness(this->header_.samples_);
		if(line.position == previous_position && *retValue == *prevRetValue){
			std::cerr << Helpers::timestamp("WARNING", "VCF") << "Duplicate position in file. Dropping... (" << this->header_[*retValue].name << ":" << line.position << ")" << std::endl;
			goto next;
		}

		// Execute only if the line is simple

		if(line.IsSimple()){
			if(missing > DEFAULT_MISSINGNESS_CUTOFF){
				std::cerr << Helpers::timestamp("WARNING", "VCF") << "Dropping " << this->header_[*retValue].name << ":" << line.position << " with " << missing*100 << "% missing values..." << std::endl;
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

	std::cerr << Helpers::timestamp("LOG", "WRITER") << "Wrote: " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.GetVariantsWritten()))
													 << " variants to " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
													 << " blocks..." << std::endl;

	// Temp
	//for(U32 i = 0; i < this->header_.getContigs(); ++i){
	//	if(this->header_.getContig(i).tomahawkBlocks > 0)
	//		std::cerr << this->header_.getContig(i) << '\t' << this->header_.getContig(i).tomahawkBlocks << std::endl;
	//}

	return true;
}

bool VCFParser::GetHeaderLines(void){
	if(!this->reader_.good())
		return false;

	if(!this->ValidateVCF())
		return false;

	// Get header lines
	while(this->reader_.getLine()){
		//std::cout << std::string(this->reader_.buffer_, this->reader_.size()) << std::endl;
		//std::cout << "last char: " << this->reader_[this->reader_.size()-2] << std::endl;
		//std::cout << "Dump size and beginning: " << this->reader_.size() << '\t' << this->reader_.buffer_[0] << '\t' << this->reader_.buffer_[1] << std::endl;

		if(this->reader_.buffer_[1] != '#')
			break;

		if(!this->header_.checkLine(this->reader_.buffer_, this->reader_.size()-2)){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Failed to validate header lines" << std::endl;
			return false;
		}

		// Todo: Parse header line

		this->reader_.clear();
	}

	std::cerr << Helpers::timestamp("LOG", "VCF") << "Parsed " << this->header_.getLines()+1 << " header lines..." << std::endl;
	return true;
}

bool VCFParser::ValidateVCF(void){
	if(!this->reader_.good())
		return false;

	if(!this->reader_.getLine()){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not validate file..." << std::endl;
		return false;
	}

	//std::cerr << std::string(&this->reader_[0], this->reader_.size()) << std::endl;
	if(strncmp(&this->reader_[0], &VCF::Constants::HEADER_VCF_FORMAT[0], VCF::Constants::HEADER_VCF_FORMAT.size()) != 0){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Invalid VCF format..." << std::endl;
		return false;
	}

	if(strncmp(&this->reader_[0], &VCF::Constants::HEADER_VCF_VERSION[0], VCF::Constants::HEADER_VCF_VERSION.size()) != 0){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Invalid VCF version < 4.x..." << std::endl;
		return false;
	}

	this->reader_.clear();

	return true;
}

bool VCFParser::SampleLine(void){
	// At broken position is main header line
	// Validate header
	if(strncmp(&Tomahawk::VCF::Constants::HEADER_COLUMN[0], &this->reader_.buffer_[0], Tomahawk::VCF::Constants::HEADER_COLUMN.size()) != 0){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not validate header line" << std::endl;
		return false;
	}

	uint32_t search_position = Tomahawk::VCF::Constants::HEADER_COLUMN.size() + 1;
	uint32_t delimiters_found = 0;
	while(true){ // while there is samples in line
		char* found = std::find(&this->reader_[search_position], &this->reader_[this->reader_.size()], Tomahawk::VCF::Constants::VCF_DELIMITER);
		if(*found != Tomahawk::VCF::Constants::VCF_DELIMITER)
			break;

		//std::cerr << std::string(&this->reader_[search_position], (found - this->reader_.buffer_ + 1) - search_position) << std::endl;
		search_position = found - this->reader_.buffer_ + 1;
		++delimiters_found;
	}
	// Last one
	//std::cerr << std::string(&this->reader_[search_position], this->reader_.size()  - search_position) << std::endl;

	std::cerr << Helpers::timestamp("LOG", "VCF") << "Found " << delimiters_found << " samples..." << std::endl;
	this->header_.setSamples(delimiters_found);

	// Parse
	search_position = Tomahawk::VCF::Constants::HEADER_COLUMN.size() + 1;
	delimiters_found = 0;
	uint32_t* retValue;
	char* found = 0;
	while(found != &this->reader_[this->reader_.size()]){ // while there is samples in line
		found = std::find(&this->reader_[search_position], &this->reader_[this->reader_.size()], Tomahawk::VCF::Constants::VCF_DELIMITER);
		//if(*found != Tomahawk::VCF::Constants::VCF_DELIMITER && found != &this->reader_[this->reader_.size()])
		//	break;

		std::string sampleName(&this->reader_[search_position], (found - this->reader_.buffer_ + 1) - search_position - 1);
		if(sampleName == "FORMAT"){
			search_position = found - this->reader_.buffer_ + 1;
			continue;
		}

		//td::cerr << sampleName << std::endl;

		if(!this->header_.getSample(sampleName, retValue))
			this->header_.addSample(sampleName);
		else {
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Duplicated sample name in header..." << std::endl;
			exit(1);
		}

		search_position = found - this->reader_.buffer_ + 1;
	}

	this->reader_.clear();
	return true;
}

}
} /* namespace Tomahawk */
