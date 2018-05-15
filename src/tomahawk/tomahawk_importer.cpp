#include <tomahawk/tomahawk_importer.h>
#include <tomahawk/tomahawk_reader.h>
#include <string>

#include "io/reader.h"
#include "import_writer.h"

namespace tomahawk {

TomahawkImporter::TomahawkImporter(const std::string inputFile, const std::string outputPrefix) :
	block_flush_limit(65536),
	input_file(inputFile),
	output_prefix(outputPrefix),
	reader_(inputFile),
	writer_(this->filters),
	vcf_header_(nullptr),
	rle_controller(nullptr)
{

}

TomahawkImporter::~TomahawkImporter(){
	delete this->rle_controller;
	// do not delete header: it might be a borrowed pointer
}

bool TomahawkImporter::Extend(std::string extendFile){

	return true;
}

bool TomahawkImporter::Build(){
	std::ifstream temp(this->input_file, std::ios::binary | std::ios::in);
	if(!temp.good()){
		std::cerr << helpers::timestamp("ERROR", "IMPORT")  << "Failed to open file..." << std::endl;
		return false;
	}
	char tempData[2];
	temp.read(&tempData[0], 2);
	temp.close();

	if(tempData[0] == '#' && tempData[1] == '#'){
		/*
		if(!this->BuildVCF()){
			std::cerr << helpers::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
			return false;
		}
		*/
		std::cerr << helpers::timestamp("ERROR", "IMPORT") << "Import file has to be binary VCF (BCF)!" << std::endl;
	} else if((BYTE)tempData[0] == io::constants::GZIP_ID1 && (BYTE)tempData[1] == io::constants::GZIP_ID2){
		if(!this->BuildBCF()){
			std::cerr << helpers::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
			return false;
		}
	} else {
		std::cerr << helpers::timestamp("ERROR", "IMPORT") << "Unknown file format!" << std::endl;
		return false;
	}
	return true;
}

bool TomahawkImporter::BuildBCF(void){
	bcf_reader_type reader;
	if(!reader.open(this->input_file)){
		std::cerr << helpers::timestamp("ERROR", "BCF")  << "Failed to open BCF file..." << std::endl;
		return false;
	}

	this->vcf_header_ = &reader.header;
	if(this->vcf_header_->samples == 0){
		std::cerr << helpers::timestamp("ERROR", "BCF") << "No samples detected in header..." << std::endl;
		return false;
	}

	if(this->vcf_header_->samples == 1){
		std::cerr << helpers::timestamp("ERROR", "IMPORT") << "Cannot run " << constants::PROGRAM_NAME << " with a single sample..." << std::endl;
		return false;
	}

	// Spawn RLE controller
	this->rle_controller = new rle_controller_type(this->vcf_header_->samples);
	this->rle_controller->DetermineBitWidth();

	this->writer_.setHeader(reader.header);
	if(!this->writer_.Open(this->output_prefix)){
		std::cerr << helpers::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	// Get a line
	bcf_entry_type entry;
	while(reader.nextVariant(entry)){
		if(entry.gt_support.hasEOV || entry.isBiallelicSimple() == false || entry.gt_support.n_missing > 3){
			std::cerr << "first: " << entry.gt_support.hasEOV << "," << entry.isBiallelicSimple() << ",miss: " << entry.gt_support.n_missing << std::endl;
			entry.reset();
			continue;
		}
		break;
	}

	if(!entry.good()){
		std::cerr << helpers::timestamp("ERROR", "IMPORT") << "No valid variants..." << std::endl;
		return false;
	}
	entry.reset();

	S32 contigID = entry.body->CHROM;
	this->sort_order_helper.previous_position  = entry.body->POS;
	this->sort_order_helper.contigID           = &contigID;
	this->sort_order_helper.prevcontigID       = contigID;
	this->writer_.totempole_entry.contigID     = contigID;
	this->writer_.totempole_entry.min_position = entry.body->POS;
	this->writer_.totempole_entry.max_position = entry.body->POS;

	if(!this->parseBCFLine(entry)){
		std::cerr << helpers::timestamp("ERROR", "BCF") << "Failed to parse BCF entry..." << std::endl;
		return false;
	}
	entry.reset();

	// Parse lines
	while(reader.nextVariant(entry)){
		if(entry.gt_support.hasEOV || entry.isBiallelicSimple() == false || entry.gt_support.n_missing > 3){
			if(entry.gt_support.hasEOV) std::cerr << "EOV: " << entry.gt_support.hasGenotypes << "," << entry.gt_support.hasEOV << "," << entry.isBiallelicSimple() << " miss: " << entry.gt_support.n_missing << std::endl;
			entry.reset();
			continue;
		}

		if(!this->parseBCFLine(entry)){
			std::cerr << helpers::timestamp("ERROR", "BCF") << "Failed to parse BCF entry..." << std::endl;
			return false;
		}

		entry.reset();
	}

	// This only happens if there are no valid entries in the file
	if(this->sort_order_helper.contigID == nullptr){
		std::cerr << helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	++this->vcf_header_->getContig(*this->sort_order_helper.contigID);
	this->writer_.flush();
	// Update container with this totempole entry
	if(this->writer_.totempole_entry.size())
		this->index += this->writer_.totempole_entry;

	this->index.buildMetaIndex(this->vcf_header_->contigs.size());
	this->writer_.WriteFinal(this->index, this->footer_);

	if(this->writer_.GetVariantsWritten() == 0){
		std::cerr << helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << helpers::timestamp("LOG", "WRITER") << "Wrote: "       << helpers::NumberThousandsSeparator(std::to_string(this->writer_.GetVariantsWritten()))
														 << " variants to " << helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
														 << " blocks..."    << std::endl;

	return true;
}

bool TomahawkImporter::BuildVCF(void){
	if(!this->reader_.open()){
		std::cerr << helpers::timestamp("ERROR","VCF") << "Failed to open file..." << std::endl;
		return false;
	}
	this->vcf_header_ = new vcf_header_type;

	if(!this->vcf_header_->parse(this->reader_)){
		std::cerr << helpers::timestamp("ERROR","VCF") << "Failed to parse VCF..." << std::endl;
		exit(1);
	}
	if(!this->vcf_header_->good()){
		std::cerr << helpers::timestamp("ERROR","VCF") << "Failed to parse VCF (" << this->vcf_header_->error_bit << ")..." << std::endl;
		return false;
	}

	if(this->vcf_header_->samples == 0){
		std::cerr << helpers::timestamp("ERROR", "IMPORT") << "No samples detected..." << std::endl;
		return false;
	}

	if(this->vcf_header_->samples == 1){
		std::cerr << helpers::timestamp("ERROR", "IMPORT") << "Cannot run " << constants::PROGRAM_NAME << " with a single sample..." << std::endl;
		return false;
	}

	// Spawn RLE controller
	this->rle_controller = new rle_controller_type(this->vcf_header_->samples);
	this->rle_controller->DetermineBitWidth();

	// Parse lines
	vcf_entry_type line(this->vcf_header_->size());

	this->reader_.clear();
	this->writer_.setHeader(*this->vcf_header_);
	if(!this->writer_.Open(this->output_prefix))
		return false;

	if(!this->reader_.getLine()){
		std::cerr << helpers::timestamp("ERROR", "VCF") << "Failed to get line" << std::endl;
		return false;
	}

	if(!line.Parse(&this->reader_[0], this->reader_.size())){
		std::cerr << helpers::timestamp("ERROR", "VCF") << "Could not parse..." << std::endl;
		return false;
	}

	// Try to get contig information from header
	if(!this->vcf_header_->getContig(std::string(line.CHROM, line.lCHROM), this->sort_order_helper.contigID)){
		std::cerr << helpers::timestamp("ERROR", "VCF") << "Contig does not exist in header..." << std::endl;
		return false;
	}
	this->sort_order_helper.prevcontigID = *this->sort_order_helper.contigID;
	this->sort_order_helper.previous_position = line.position;
	this->writer_.totempole_entry.contigID = *this->sort_order_helper.contigID;

	if(!this->parseVCFLine(line)){
		std::cerr << helpers::timestamp("ERROR", "VCF") << "Failed parse" << std::endl;
		return false;
	}

	// While there are lines
	while(this->reader_.getLine()){
		// Parse them
		if(!this->parseVCFLine(line)){
			return false;
		}
	} // end while there are vcf lines

	// This only happens if there are no valid entries in the file
	if(this->sort_order_helper.contigID == nullptr){
		std::cerr << helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	++this->vcf_header_->getContig(*this->sort_order_helper.contigID);
	this->writer_.flush();
	// Update container with this totempole entry
	if(this->writer_.totempole_entry.size())
		this->index += this->writer_.totempole_entry;

	//		return false;
	this->index.buildMetaIndex(this->vcf_header_->contigs.size());
	this->writer_.WriteFinal(this->index, this->footer_);

	if(this->writer_.GetVariantsWritten() == 0){
		std::cerr << helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << helpers::timestamp("LOG", "WRITER") << "Wrote: " << helpers::NumberThousandsSeparator(std::to_string(this->writer_.GetVariantsWritten()))
														 << " variants to " << helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
														 << " blocks..." << std::endl;

	delete this->vcf_header_;

	return true;
}

bool TomahawkImporter::parseBCFLine(bcf_entry_type& line){
	if(this->sort_order_helper.prevcontigID != line.body->CHROM){
		if(line.body->CHROM < this->sort_order_helper.prevcontigID){
			std::cerr << helpers::timestamp("ERROR", "IMPORT") << "Contigs are not sorted (" << (*this->vcf_header_)[this->sort_order_helper.prevcontigID].name << " > " << (*this->vcf_header_)[line.body->CHROM].name << ")..." << std::endl;
			exit(1);
		}

		if(!SILENT){
			std::cerr << this->writer_.totempole_entry.n_variants << std::endl;
			std::cerr << helpers::timestamp("LOG", "IMPORT") << "Switch detected: " << this->vcf_header_->getContig(this->sort_order_helper.prevcontigID).name << "->" << this->vcf_header_->getContig(line.body->CHROM).name << "..." << std::endl;
		}

		this->sort_order_helper.previous_position = 0;

		// Get new contig value from header
		// and flush out data
		++this->vcf_header_->getContig(line.body->CHROM);
		this->writer_.flush();

		// Update container with this totempole entry
		if(this->writer_.totempole_entry.size())
			this->index += this->writer_.totempole_entry;

		// Update index values
		this->writer_.TotempoleSwitch(line.body->CHROM, 0);
	}

	// Assert position is in range
	if(line.body->POS + 1 > this->vcf_header_->getContig(line.body->CHROM).length){
		std::cerr << helpers::timestamp("ERROR", "IMPORT") << (*this->vcf_header_)[line.body->CHROM].name << ':' << line.body->POS+1 << " > reported max size of contig (" << (*this->vcf_header_)[line.body->CHROM].length << ")..." << std::endl;
		return false;
	}

	// Assert file is ordered
	if(line.body->POS + 1 < this->sort_order_helper.previous_position){
		std::cerr << helpers::timestamp("ERROR", "IMPORT") << "File is not sorted by coordinates (" << (*this->vcf_header_)[line.body->CHROM].name << ':' << line.body->POS+1 << " > " << (*this->vcf_header_)[line.body->CHROM].name << ':' << this->sort_order_helper.previous_position << ")..." << std::endl;
		return false;
	}


	// Assess missingness
	if(line.body->POS == this->sort_order_helper.previous_position && line.body->CHROM == this->sort_order_helper.prevcontigID){
		if(this->sort_order_helper.previous_included){
			//if(!SILENT)
			//	std::cerr << helpers::timestamp("WARNING", "BCF") << "Duplicate position (" << (*this->header_)[line.body->CHROM].name << ":" << line.body->POS+1 << "): Dropping..." << std::endl;

			goto next;
		} else {
			//if(!SILENT)
			//	std::cerr << helpers::timestamp("WARNING", "BCF") << "Duplicate position (" << (*this->header_)[line.body->CHROM].name << ":" << line.body->POS+1 << "): Keeping (drop other)..." << std::endl;

		}
	}

	// Execute only if the line is simple (biallelic and SNP)
	if(line.isSimple()){
		// Flush if output block is over some size
		if(this->writer_.checkSize()){
			++this->vcf_header_->getContig(line.body->CHROM); // update block count for this contigID
			this->writer_.flush();

			// Update container with this totempole entry
			if(this->writer_.totempole_entry.size())
				this->index += this->writer_.totempole_entry;

			this->writer_.TotempoleSwitch(line.body->CHROM, 0);
		}
		if(this->writer_.add(line))
			this->sort_order_helper.previous_included = true;
		else
			this->sort_order_helper.previous_included = false;
	} else
		this->sort_order_helper.previous_included = false;

	next:
	this->sort_order_helper.previous_position = line.body->POS + 1;
	this->sort_order_helper.prevcontigID      = line.body->CHROM;

	return true;
}

bool TomahawkImporter::parseVCFLine(vcf_entry_type& line){
	// Parse a VCF line
	if(!line.Parse(&this->reader_[0], this->reader_.size())){
		std::cerr << helpers::timestamp("ERROR", "VCF") << "Could not parse..." << std::endl;
		return false;
	}

	// Try to get contig information from header
	if(!this->vcf_header_->getContig(std::string(line.CHROM, line.lCHROM), this->sort_order_helper.contigID)){
		std::cerr << helpers::timestamp("ERROR", "VCF") << "Contig does not exist in header..." << std::endl;
		return false;
	}

	// Switch in chromosome detected
	if(this->sort_order_helper.prevcontigID != *this->sort_order_helper.contigID){
		if(*this->sort_order_helper.contigID < this->sort_order_helper.prevcontigID){
			std::cerr << helpers::timestamp("ERROR", "VCF") << "Contigs are not sorted (" << (*this->vcf_header_)[this->sort_order_helper.prevcontigID].name << " > " << (*this->vcf_header_)[*this->sort_order_helper.contigID].name << ")..." << std::endl;
			exit(1);
		}

		if(!SILENT)
			std::cerr << helpers::timestamp("LOG", "VCF") << "Switch detected: " << this->vcf_header_->getContig(this->sort_order_helper.prevcontigID).name << "->" << this->vcf_header_->getContig(*this->sort_order_helper.contigID).name << "..." << std::endl;

		this->sort_order_helper.previous_position = 0;

		// Get new contig value from header
		// and flush out data
		++this->vcf_header_->getContig(*this->sort_order_helper.contigID);
		this->writer_.flush();

		// Update container with this totempole entry
		if(this->writer_.totempole_entry.size())
			this->index += this->writer_.totempole_entry;

		// Update index values
		this->writer_.TotempoleSwitch(*this->sort_order_helper.contigID, 0);
	}

	// Assert position is in range
	if(line.position > this->vcf_header_->getContig(*this->sort_order_helper.contigID).length){
		std::cerr << helpers::timestamp("ERROR", "VCF") << (*this->vcf_header_)[*this->sort_order_helper.contigID].name << ':' << line.position << " > reported max size of contig (" << (*this->vcf_header_)[*this->sort_order_helper.contigID].length << ")..." << std::endl;
		return false;
	}

	// Assert file is ordered
	if(line.position < this->sort_order_helper.previous_position){
		std::cerr << helpers::timestamp("ERROR", "VCF") << "File is not sorted by coordinates (" << (*this->vcf_header_)[*this->sort_order_helper.contigID].name << ':' << line.position << " > " << (*this->vcf_header_)[*this->sort_order_helper.contigID].name << ':' << this->sort_order_helper.previous_position << ")..." << std::endl;
		return false;
	}

	// Execute only if the line is simple (biallelic and SNP)
	if(line.IsSimple()){
		// Only check missing if simple
		const double missing = line.getMissingness(this->vcf_header_->samples);
		if(line.position == this->sort_order_helper.previous_position && *this->sort_order_helper.contigID == this->sort_order_helper.prevcontigID){
			if(this->sort_order_helper.previous_included){
				//if(!SILENT)
				//	std::cerr << helpers::timestamp("WARNING", "VCF") << "Duplicate position (" << (*this->header_)[*this->sort_order_helper.contigID].name << ":" << line.position << "): Dropping..." << std::endl;

				goto next;
			} else {
				//if(!SILENT)
				//	std::cerr << helpers::timestamp("WARNING", "VCF") << "Duplicate position (" << (*this->header_)[*this->sort_order_helper.contigID].name << ":" << line.position << "): Keeping (drop other)..." << std::endl;
			}
		}

		if(missing > this->filters.missingness){
			//if(!SILENT)
			//	std::cerr << helpers::timestamp("WARNING", "VCF") << "Large missingness (" << (*this->header_)[*this->sort_order_helper.contigID].name << ":" << line.position << ", " << missing*100 << "%).  Dropping..." << std::endl;

			this->sort_order_helper.previous_included = false;
			goto next;
		}

		// Flush if output block is over some size
		if(this->writer_.checkSize()){
			++this->vcf_header_->getContig(*this->sort_order_helper.contigID); // update block count for this contigID
			this->writer_.flush();

			// Update container with this totempole entry
			if(this->writer_.totempole_entry.size())
				this->index += this->writer_.totempole_entry;

			this->writer_.TotempoleSwitch(*this->sort_order_helper.contigID, this->sort_order_helper.previous_position);
		}
		if(this->writer_.add(line))
			this->sort_order_helper.previous_included = true;
		else
			this->sort_order_helper.previous_included = false;
	} else
		this->sort_order_helper.previous_included = false;

	next:
	this->sort_order_helper.previous_position = line.position;
	this->sort_order_helper.prevcontigID = *this->sort_order_helper.contigID;
	this->reader_.clear();
	return true;
}

} /* namespace Tomahawk */
