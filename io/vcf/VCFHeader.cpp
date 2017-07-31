#include "VCFHeader.h"

namespace Tomahawk {
namespace VCF{

VCFHeader::VCFHeader() :
	error_bit(VCF_PASS),
	samples(0),
	version(0),
	contigsHashTable(nullptr),
	sampleHashTable(nullptr)
{
}

VCFHeader::~VCFHeader(){
	delete this->contigsHashTable;
	delete this->sampleHashTable;
}

bool VCFHeader::parse(reader& stream){
	if(!this->__parseFirstLine(stream))
		return false;

	// Read remainder lines
	if(!this->__parseHeaderLines(stream))
		return false;

	if(!this->buildContigTable())
		return false;

	// Read samples line
	if(!this->__parseSampleLine(stream))
		return false;

	return true;
}

void VCFHeader::buildSampleTable(U64 samples){
	this->samples = samples;
	if(this->sampleHashTable != nullptr)
		delete this->sampleHashTable;

	if(this->samples < 1024)
		this->sampleHashTable = new hash_table(1024);
	else
		this->sampleHashTable = new hash_table(this->samples * 2);
}

bool VCFHeader::checkLine(const char* data, const U32 length){
	VCFHeaderLine line(data, length);
	if(line.Parse()){
		//std::cerr << "Input line: " << line << std::endl;

		// If the line is a contig line: make sure it is legal
		// for our purposes
		if(line.isCONTIG()){
			contig_type contig;
			BYTE found = 0;

			// Contig line has two values:
			// ID: contig name
			// length: for length is bp
			for(U32 i = 0; i < line.size(); ++i){
				if(strncmp(&line[i].KEY[0], "ID", 2) == 0){
					contig.name = std::string(line[i].VALUE);
					//std::cerr << std::string(line[i].VALUE, line[i].lVALUE) << std::endl;
					++found;
				} else if(strncmp(&line[i].KEY[0], "length", 6) == 0){
					contig.length = atoi(&line[i].VALUE[0]);
					//std::cerr << contig.length << std::endl;
					++found;
				}
			}

			// Throw error if this pattern is not found
			if(found != 2){
				std::cerr << Helpers::timestamp("WARNING","VCF") << "Illegal contig entry line with no length defined!" << std::endl;
				std::cerr << Helpers::timestamp("WARNING","VCF") << "Offending line: " << std::string(data, length+1) << std::endl;
				contig.length = std::numeric_limits<U32>::max();
			}
			this->contigs.push_back(contig);
		}

		this->lines.push_back(line); // parseable lines
		this->literal_lines.push_back(std::string(data, length + 1));
		return true;
	}

	std::cerr << Helpers::timestamp("ERROR","VCF") << "Failed to parse VCF..." << std::endl;
	return false;
}

bool VCFHeader::buildContigTable(void){
	U32* retValue;

	if(this->contigsHashTable)
		delete this->contigsHashTable;

	if(this->contigs.size() < 1024)
		this->contigsHashTable = new hash_table(1024);
	else
		this->contigsHashTable = new hash_table(this->contigs.size() * 2);

	std::cerr << Helpers::timestamp("LOG", "VCF") << "Constructing lookup table for " << this->contigs.size() << " contigs..." << std::endl;

	for(U32 i = 0; i < this->contigs.size(); ++i){
		//std::cerr << this->contigs[i] << std::endl;
		if(!(*this).getContig(this->contigs[i].name, retValue)){
			(*this).addContig(this->contigs[i].name, i);
		} else {
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Duplicated contig found (" << this->getContig(*retValue).name << "). Illegal..." << std::endl;
			this->error_bit = VCF_ERROR_LINES;
			return false;
		}
	}
	return true;
}

bool VCFHeader::__parseFirstLine(reader& stream){
	if(!stream.good()){
		this->error_bit = STREAM_BAD;
		return false;
	}

	if(!stream.getLine()){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not validate file..." << std::endl;
		this->error_bit = STREAM_BAD;
		return false;
	}

	// Parse
	if(strncmp(&stream[0], &VCF::Constants::HEADER_VCF_FORMAT[0], VCF::Constants::HEADER_VCF_FORMAT.size()) != 0){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Invalid VCF format..." << std::endl;
		this->error_bit = VCF_ERROR_LINE1;
		return false;
	}

	if(strncmp(&stream[0], &VCF::Constants::HEADER_VCF_VERSION[0], VCF::Constants::HEADER_VCF_VERSION.size()) != 0){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Invalid VCF version < 4.x..." << std::endl;
		this->error_bit = VCF_ERROR_LINE1;
		return false;
	}

	stream.clear();
	return true;
}

bool VCFHeader::__parseHeaderLines(reader& stream){
	while(stream.getLine()){
		if(stream.buffer_[1] != '#')
			break;

		if(!this->checkLine(stream.buffer_, stream.size() - 2)){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Failed to validate header lines" << std::endl;
			this->error_bit = VCF_ERROR_LINES;
			return false;
		}

		stream.clear();
	}

	if(!stream.good()){
		this->error_bit = STREAM_BAD;
		return false;
	}

	return true;
}

bool VCFHeader::__parseSampleLine(reader& stream){
	// At broken position is main header line
	// Validate header
	if(strncmp(&Tomahawk::VCF::Constants::HEADER_COLUMN[0], &stream.buffer_[0], Tomahawk::VCF::Constants::HEADER_COLUMN.size()) != 0){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not validate header line" << std::endl;
		this->error_bit = VCF_ERROR_SAMPLE;
		return false;
	}

	U32 search_position = Tomahawk::VCF::Constants::HEADER_COLUMN.size() + 1;
	U64 delimiters_found = 0;
	while(true){ // while there is samples in line
		char* found = std::find(&stream[search_position], &stream[stream.size()], Tomahawk::VCF::Constants::VCF_DELIMITER);
		if(*found != Tomahawk::VCF::Constants::VCF_DELIMITER)
			break;

		//std::cerr << std::string(&stream[search_position], (found - stream.buffer_ + 1) - search_position) << std::endl;
		search_position = found - stream.buffer_ + 1;
		++delimiters_found;
	}

	this->buildSampleTable(delimiters_found);

	// Parse
	search_position = Tomahawk::VCF::Constants::HEADER_COLUMN.size() + 1;
	delimiters_found = 0;
	U32* retValue;
	char* found = 0;
	while(found != &stream[stream.size()]){ // while there are samples in line
		found = std::find(&stream[search_position], &stream[stream.size()], Tomahawk::VCF::Constants::VCF_DELIMITER);
		//if(*found != Tomahawk::VCF::Constants::VCF_DELIMITER && found != &stream[stream.size()])
		//	break;

		std::string sampleName(&stream[search_position], (found - stream.buffer_ + 1) - search_position - 1);
		if(sampleName == "FORMAT"){
			search_position = found - stream.buffer_ + 1;
			continue;
		}

		//std::cerr << sampleName << std::endl;

		if(!this->getSample(sampleName, retValue))
			this->addSample(sampleName);
		else {
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Duplicated sample name in header..." << std::endl;
			this->error_bit = VCF_ERROR_LINES;
		}

		search_position = found - stream.buffer_ + 1;
	}

	stream.clear();
	return true;
}


}
} /* namespace Tomahawk */
