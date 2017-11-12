#include "TomahawkImporter.h"
#include "TomahawkImportWriter.h"
#include "TomahawkReader.h"

namespace Tomahawk {

TomahawkImporter::TomahawkImporter(std::string inputFile, std::string outputPrefix, const U32 checkpoint) :
	checkpoint_size(checkpoint),
	block_flush_limit(65536),
	inputFile(inputFile),
	outputPrefix(outputPrefix),
	reader_(inputFile),
	writer_(this->filters),
	header_(nullptr),
	encoder(nullptr)
{}

TomahawkImporter::~TomahawkImporter(){
	delete this->encoder;
	// do not delete header: it might be a borrowed pointer
}

bool TomahawkImporter::Extend(std::string extendFile){
	if(this->inputFile.size() == 0){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "No input file provided..." << std::endl;
		return false;
	}

	if(extendFile.size() == 0){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "No file to extend provided..." << std::endl;
		return false;
	}

	if(!this->reader_.open()){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Failed to open file..." << std::endl;
		return false;
	}

	TomahawkReader tReader;
	if(!tReader.Open(extendFile)){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") <<  "Failed to read file..." << std::endl;
		return false;
	}

	const Totempole::TotempoleReader& totempole = tReader.getTotempole();
	*this->header_ = totempole; // Convert data in totempole to VCF header

	// Parse lines
	line_type line(totempole.getHeader().samples);

	// Spawn RLE controller
	this->encoder = new encoder_type(this->header_->samples);

	this->reader_.clear();
	// seek reader until line does not start with '#'
	std::string templine;
	while(getline(this->reader_.stream_, templine)){
		if(templine[0] != '#')
			break;
	}
	this->reader_.stream_.seekg((U64)this->reader_.stream_.tellg() - templine.size() - 1);

	this->sort_order_helper.previous_position = totempole.back().maxPosition;
	this->sort_order_helper.prevcontigID = totempole.back().contigID;

	this->writer_.setHeader(*this->header_);
	this->writer_.n_blocksWritten = totempole.getHeader().blocks;
	this->writer_.largest_uncompressed_block = totempole.getHeader().largest_uncompressed;
	if(!this->writer_.OpenExtend(extendFile))
		return false;

	// While there are lines
	while(this->reader_.getLine()){
		// Parse them
		if(!this->parseVCFLine(line)){
			return false;
		}
	} // end while there are vcf lines

	// This only happens if there are no valid entries in the file
	if(this->sort_order_helper.contigID == nullptr){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	++this->header_->getContig(*this->sort_order_helper.contigID);
	this->writer_.flush(this->permutator);

	this->writer_.WriteFinal();

	if(this->writer_.getVariantsWritten() == 0){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Wrote: " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.getVariantsWritten()))
														 << " variants to " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
														 << " blocks..." << std::endl;

	// Garbage
	this->header_->unsetBorrowedPointers();
	return true;
}

bool TomahawkImporter::Build(){
	std::ifstream temp(this->inputFile, std::ios::binary | std::ios::in);
	if(!temp.good()){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT")  << "Failed to open file..." << std::endl;
		return false;
	}
	char tempData[2];
	temp.read(&tempData[0], 2);
	temp.close();

	if(tempData[0] == '#' && tempData[1] == '#'){
		if(!this->BuildVCF()){
			std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
			return false;
		}
	} else if((BYTE)tempData[0] == IO::Constants::GZIP_ID1 && (BYTE)tempData[1] == IO::Constants::GZIP_ID2){
		if(!this->BuildBCF()){
			std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
			return false;
		}
	} else {
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Unknown file format!" << std::endl;
		return false;
	}
	return true;
}

bool TomahawkImporter::BuildBCF(void){
	bcf_reader_type reader;
	if(!reader.open(this->inputFile)){
		std::cerr << Helpers::timestamp("ERROR", "BCF")  << "Failed to open BCF file..." << std::endl;
		return false;
	}

	this->header_ = &reader.header;
	if(this->header_->samples == 0){
		std::cerr << Helpers::timestamp("ERROR", "BCF") << "No samples detected in header..." << std::endl;
		return false;
	}

	if(this->header_->samples == 1){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Cannot run " << Tomahawk::Constants::PROGRAM_NAME << " with a single sample..." << std::endl;
		return false;
	}

	// Spawn RLE controller
	this->encoder = new encoder_type(this->header_->samples);
	//this->encoder->DetermineBitWidth();
	this->permutator.setSamples(this->header_->samples);

	this->writer_.setHeader(reader.header);
	if(!this->writer_.Open(this->outputPrefix)){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	///
	/// TODO
	/// temp

	this->sort_order_helper.previous_position = 0;
	this->sort_order_helper.contigID = nullptr;
	this->sort_order_helper.prevcontigID = 0;
	this->writer_.totempole_entry.contigID = 0;
	this->writer_.totempole_entry.minPosition = 0;

	//std::cerr << "PPA_conventional\tPPA_best\tPPA_byte\tPPA_u16\tPPA_u32\tPPA_u64\trle_conventional\trle_best\trle_byte\trle_u16\trle_u32\trle_u64\tfd_rle_best_ppa_best\tmemory_savings_rle_ppa\tfc_rle_conventional_ppa_best" << std::endl;
	while(true){
		if(!reader.getVariants(this->checkpoint_size))
			break;

		S32 contigID = reader[0].body->CHROM;
		this->sort_order_helper.previous_position = reader[0].body->POS;
		this->sort_order_helper.contigID = &contigID;
		this->sort_order_helper.prevcontigID = contigID;
		this->writer_.totempole_entry.contigID = contigID;
		this->writer_.totempole_entry.minPosition = reader[0].body->POS;

		// Reset permutate
		if(!this->permutator.build(reader)){
			std::cerr << "fail" << std::endl;
			return false;
		}

		for(U32 i = 0; i < reader.size(); ++i){
			if(!this->parseBCFLine(reader[i])){
				std::cerr << "failed to parse" << std::endl;
				return false;
			}
		}

		++this->header_->getContig(contigID); // update block count for this contigID
		this->writer_.flush(this->permutator);
		this->writer_.TotempoleSwitch(contigID, this->sort_order_helper.previous_position);

		// Reset permutator
		this->permutator.reset();
	}

	// This only happens if there are no valid entries in the file
	if(this->sort_order_helper.contigID == nullptr){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	++this->header_->getContig(*this->sort_order_helper.contigID);
	this->writer_.flush(this->permutator);
	this->writer_.WriteFinal();

	if(this->writer_.getVariantsWritten() == 0){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Wrote: " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.getVariantsWritten()))
														 << " variants to " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
														 << " blocks..." << std::endl;

	return(true);
}

bool TomahawkImporter::BuildVCF(void){
	if(!this->reader_.open()){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Failed to open file..." << std::endl;
		return false;
	}
	this->header_ = new header_type;

	if(!this->header_->parse(this->reader_)){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Failed to parse VCF..." << std::endl;
		exit(1);
	}
	if(!this->header_->good()){
		std::cerr << Helpers::timestamp("ERROR","VCF") << "Failed to parse VCF (" << this->header_->error_bit << ")..." << std::endl;
		return false;
	}

	if(this->header_->samples == 0){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "No samples detected..." << std::endl;
		return false;
	}

	if(this->header_->samples == 1){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Cannot run " << Tomahawk::Constants::PROGRAM_NAME << " with a single sample..." << std::endl;
		return false;
	}

	// Spawn RLE controller
	this->encoder = new encoder_type(this->header_->samples);
	//this->encoder->DetermineBitWidth();

	// Parse lines
	line_type line(this->header_->size());

	this->reader_.clear();
	this->writer_.setHeader(*this->header_);
	if(!this->writer_.Open(this->outputPrefix))
		return false;

	if(!this->reader_.getLine()){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Failed to get line" << std::endl;
		return false;
	}

	if(!line.Parse(&this->reader_[0], this->reader_.size())){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not parse..." << std::endl;
		return false;
	}

	// Try to get contig information from header
	if(!this->header_->getContig(std::string(line.CHROM, line.lCHROM), this->sort_order_helper.contigID)){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Contig does not exist in header..." << std::endl;
		return false;
	}
	this->sort_order_helper.prevcontigID = *this->sort_order_helper.contigID;
	this->sort_order_helper.previous_position = line.position;
	this->writer_.totempole_entry.contigID = *this->sort_order_helper.contigID;

	if(!this->parseVCFLine(line)){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Failed parse" << std::endl;
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
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	++this->header_->getContig(*this->sort_order_helper.contigID);
	this->writer_.flush(this->permutator);
	//		return false;

	this->writer_.WriteFinal();

	if(this->writer_.getVariantsWritten() == 0){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Wrote: " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.getVariantsWritten()))
														 << " variants to " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
														 << " blocks..." << std::endl;

	delete this->header_;

	return true;
}

bool TomahawkImporter::parseBCFLine(bcf_entry_type& line){
	if(this->sort_order_helper.prevcontigID != line.body->CHROM){
		if(line.body->CHROM < this->sort_order_helper.prevcontigID){
			std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Contigs are not sorted (" << (*this->header_)[this->sort_order_helper.prevcontigID].name << " > " << (*this->header_)[line.body->CHROM].name << ")..." << std::endl;
			exit(1);
		}

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "IMPORT") << "Switch detected: " << this->header_->getContig(this->sort_order_helper.prevcontigID).name << "->" << this->header_->getContig(line.body->CHROM).name << "..." << std::endl;

		this->sort_order_helper.previous_position = 0;

		// Get new contig value from header
		// and flush out data
		++this->header_->getContig(line.body->CHROM);
		//this->writer_.flush(this->permutator);

		// Update index values
		this->writer_.TotempoleSwitch(line.body->CHROM, 0);
	}

	// Assert position is in range
	if(line.body->POS + 1 > this->header_->getContig(line.body->CHROM).length){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << (*this->header_)[line.body->CHROM].name << ':' << line.body->POS+1 << " > reported max size of contig (" << (*this->header_)[line.body->CHROM].length << ")..." << std::endl;
		return false;
	}

	// Assert file is ordered
	if(line.body->POS < this->sort_order_helper.previous_position){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "File is not sorted by coordinates (" << (*this->header_)[line.body->CHROM].name << ':' << line.body->POS+1 << " > " << (*this->header_)[line.body->CHROM].name << ':' << this->sort_order_helper.previous_position << ")..." << std::endl;
		return false;
	}

	// Assess missingness
	const double missing = line.getMissingness(this->header_->samples);

	// Execute only if the line is simple (biallelic and SNP)
	if(line.isSimple()){
		if(missing > this->filters.missingness){
			//if(!SILENT)
			//std::cerr << Helpers::timestamp("WARNING", "VCF") << "Large missingness (" << (*this->header_)[line.body->CHROM].name << ":" << line.body->POS+1 << ", " << missing << "%).  Dropping... / " << this->filters.missingness << std::endl;
			goto next;
		}

		// Flush if output block is over some size
		/*
		if(this->writer_.checkSize()){
			++this->header_->getContig(line.body->CHROM); // update block count for this contigID
			this->writer_.flush();

			this->writer_.TotempoleSwitch(line.body->CHROM, this->sort_order_helper.previous_position);
		}
		*/
		this->writer_.add(line, this->permutator.getPPA());
	}

	next:
	this->sort_order_helper.previous_position = line.body->POS;
	this->sort_order_helper.prevcontigID = line.body->CHROM;

	return true;
}

bool TomahawkImporter::parseVCFLine(line_type& line){
	// Parse a VCF line
	if(!line.Parse(&this->reader_[0], this->reader_.size())){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not parse..." << std::endl;
		return false;
	}

	// Try to get contig information from header
	if(!this->header_->getContig(std::string(line.CHROM, line.lCHROM), this->sort_order_helper.contigID)){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Contig does not exist in header..." << std::endl;
		return false;
	}

	// Switch in chromosome detected
	if(this->sort_order_helper.prevcontigID != *this->sort_order_helper.contigID){
		if(*this->sort_order_helper.contigID < this->sort_order_helper.prevcontigID){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Contigs are not sorted (" << (*this->header_)[this->sort_order_helper.prevcontigID].name << " > " << (*this->header_)[*this->sort_order_helper.contigID].name << ")..." << std::endl;
			exit(1);
		}

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "VCF") << "Switch detected: " << this->header_->getContig(this->sort_order_helper.prevcontigID).name << "->" << this->header_->getContig(*this->sort_order_helper.contigID).name << "..." << std::endl;

		this->sort_order_helper.previous_position = 0;

		// Get new contig value from header
		// and flush out data
		++this->header_->getContig(*this->sort_order_helper.contigID);
		this->writer_.flush(this->permutator);

		// Update index values
		this->writer_.TotempoleSwitch(*this->sort_order_helper.contigID, 0);
	}

	// Assert position is in range
	if(line.position > this->header_->getContig(*this->sort_order_helper.contigID).length){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << (*this->header_)[*this->sort_order_helper.contigID].name << ':' << line.position << " > reported max size of contig (" << (*this->header_)[*this->sort_order_helper.contigID].length << ")..." << std::endl;
		return false;
	}

	// Assert file is ordered
	if(line.position < this->sort_order_helper.previous_position){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "File is not sorted by coordinates (" << (*this->header_)[*this->sort_order_helper.contigID].name << ':' << line.position << " > " << (*this->header_)[*this->sort_order_helper.contigID].name << ':' << this->sort_order_helper.previous_position << ")..." << std::endl;
		return false;
	}

	// Execute only if the line is simple (biallelic and SNP)
	if(line.IsSimple()){
		// Only check missing if simple
		const double missing = line.getMissingness(this->header_->samples);

		if(missing > this->filters.missingness){
			//if(!SILENT)
			//	std::cerr << Helpers::timestamp("WARNING", "VCF") << "Large missingness (" << (*this->header_)[*this->sort_order_helper.contigID].name << ":" << line.position << ", " << missing*100 << "%).  Dropping..." << std::endl;
			goto next;
		}

		// Flush if output block is over some size
		if(this->writer_.checkSize()){
			++this->header_->getContig(*this->sort_order_helper.contigID); // update block count for this contigID
			this->writer_.flush(this->permutator);

			this->writer_.TotempoleSwitch(*this->sort_order_helper.contigID, this->sort_order_helper.previous_position);
		}
		this->writer_.add(line);
	}

	next:
	this->sort_order_helper.previous_position = line.position;
	this->sort_order_helper.prevcontigID = *this->sort_order_helper.contigID;
	this->reader_.clear();
	return true;
}

} /* namespace Tomahawk */
