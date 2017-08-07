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

VCFParser::~VCFParser(){
	delete this->rle_controller;
}

bool VCFParser::Extend(){
	// reader open
	// header read
	// header validate
	// twi open
	// twi read
	// twi header validate
	// seek eof
	// check eof
	// read last block
	// get last contigID and last position
	return false;
}

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

	// Spawn RLE controller
	this->rle_controller = new rle_controller_type(this->header_.samples);
	this->rle_controller->DetermineBitWidth();

	// Parse lines
	line_type line(this->header_.size());

	this->reader_.clear();
	this->writer_.setHeader(this->header_);
	// Todo
	if(!this->writer_.Open(this->outputPrefix))
		return false;

	if(!this->reader_.getLine()){
		std::cerr << "failed to get line" << std::endl;
		return false;
	}
	// Parse a VCF line
	if(!line.Parse(&this->reader_[0], this->reader_.size())){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not parse..." << std::endl;
		return false;
	}

	this->sort_order_helper.previous_position = line.position;
	// Try to get contig information from header
	if(!this->header_.getContig(std::string(line.CHROM, line.lCHROM), this->sort_order_helper.contigID)){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Contig does not exist in header..." << std::endl;
		return false;
	}
	this->sort_order_helper.prevcontigID = this->sort_order_helper.contigID;
	if(!this->parseLine(line)){
		std::cerr << "faiaeld parse" << std::endl;
		return false;
	}

	// While there are lines
	while(this->reader_.getLine()){
		// Parse them
		if(!this->parseLine(line)){
			return false;
		}
	} // end while there are vcf lines

	// This only happens if there are no valid entries in the file
	if(this->sort_order_helper.contigID == nullptr){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Did not import any variants..." << std::endl;
		return false;
	}

	++this->header_.getContig(*this->sort_order_helper.contigID);
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

bool VCFParser::parseLine(line_type& line){
	// Parse a VCF line
	if(!line.Parse(&this->reader_[0], this->reader_.size())){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not parse..." << std::endl;
		return false;
	}

	// Try to get contig information from header
	if(!this->header_.getContig(std::string(line.CHROM, line.lCHROM), this->sort_order_helper.contigID)){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Contig does not exist in header..." << std::endl;
		return false;
	}

	// Switch in chromosome detected
	if(this->sort_order_helper.prevcontigID != this->sort_order_helper.contigID){
		if(*this->sort_order_helper.contigID < *this->sort_order_helper.prevcontigID){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Contigs are not sorted (" << this->header_[*this->sort_order_helper.prevcontigID].name << " > " << this->header_[*this->sort_order_helper.contigID].name << ")..." << std::endl;
			exit(1);
		}

		//std::cerr << Helpers::timestamp("DEBUG", "VCF") << "Switch detected: " << this->header_.getContig(*prevcontigID).name << "->" << this->header_.getContig(*contigID).name << "..." << std::endl;
		this->sort_order_helper.previous_position = 0;

		// Get new contig value from header
		// and flush out data
		++this->header_.getContig(*this->sort_order_helper.contigID);
		this->writer_.flush();

		// Update index values
		this->writer_.TotempoleSwitch(*this->sort_order_helper.contigID, 0);
	}

	// Assert position is in range
	if(line.position > this->header_.getContig(*this->sort_order_helper.contigID).length){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << this->header_[*this->sort_order_helper.contigID].name << ':' << line.position << " > reported max size of contig (" << this->header_[*this->sort_order_helper.contigID].length << ")..." << std::endl;
		return false;
	}

	// Assert file is ordered
	if(line.position < this->sort_order_helper.previous_position){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "File is not sorted by coordinates (" << this->header_[*this->sort_order_helper.contigID].name << ':' << line.position << " > " << this->header_[*this->sort_order_helper.contigID].name << ':' << this->sort_order_helper.previous_position << ")..." << std::endl;
		return false;
	}

	// Assess missingness
	const float missing = line.getMissingness(this->header_.samples);
	if(line.position == this->sort_order_helper.previous_position && this->sort_order_helper.contigID == this->sort_order_helper.prevcontigID){
		if(!SILENT)
			std::cerr << Helpers::timestamp("WARNING", "VCF") << "Duplicate position (" << this->header_[*this->sort_order_helper.contigID].name << ":" << line.position << "): Dropping..." << std::endl;

		goto next;
	}

	// Execute only if the line is simple (biallelic and SNP)
	if(line.IsSimple()){
		if(missing > DEFAULT_MISSINGNESS_CUTOFF){
			if(!SILENT)
				std::cerr << Helpers::timestamp("WARNING", "VCF") << "Large missingness (" << this->header_[*this->sort_order_helper.contigID].name << ":" << line.position << ", " << missing*100 << "%).  Dropping..." << std::endl;

			goto next;
		}

		// Flush if output block is over some size
		if(this->writer_.checkSize()){
			++this->header_.getContig(*this->sort_order_helper.contigID); // update block count for this contigID
			this->writer_.flush();

			this->writer_.TotempoleSwitch(*this->sort_order_helper.contigID, this->sort_order_helper.previous_position);
		}
		this->writer_ += line;
	}

	next:
	this->sort_order_helper.previous_position = line.position;
	this->sort_order_helper.prevcontigID = this->sort_order_helper.contigID;
	this->reader_.clear();
	return true;
}

}
} /* namespace Tomahawk */
