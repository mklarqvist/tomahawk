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

	// Copy string literal to header
	U32 curPos = stream.tellg();
	U32 headerLength = curPos - stream.size();
	this->literal.resize(headerLength);
	stream.stream_.seekg(0);
	stream.stream_.read(&this->literal[0], headerLength);
	stream.stream_.seekg(curPos);

	// Read samples line
	if(!this->__parseSampleLine(stream))
		return false;

	return true;
}

bool VCFHeader::parse(const char* const data, const U32& length){
	U32 offset = 0;
	if(!this->__parseFirstLine(data, offset))
		return false;

	// Read remainder lines
	if(!this->__parseHeaderLines(data, offset))
		return false;

	if(!this->buildContigTable())
		return false;

	// Read samples line
	if(!this->__parseSampleLine(data, offset, length))
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
		// If the line is a contig line: make sure it is legal
		// for our purposes
		if(line.isCONTIG()){
			contig_type contig;
			BYTE found = 0;

			// Contig line has two values:
			// ID: contig name
			// length: for length is bp
			for(U32 i = 0; i < line.size(); ++i){
				if(strncmp(&line[i].KEY[0], "ID", 2) == 0 && line[i].KEY.size() == 2){
					contig.name = std::string(line[i].VALUE);
					//std::cerr << std::string(line[i].VALUE, line[i].lVALUE) << std::endl;
					++found;
				} else if(strncmp(&line[i].KEY[0], "length", 6) == 0 && line[i].KEY.size() == 6){
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
	S32* retValue;

	if(this->contigsHashTable)
		delete this->contigsHashTable;

	if(this->contigs.size() < 1024)
		this->contigsHashTable = new hash_table(1024);
	else
		this->contigsHashTable = new hash_table(this->contigs.size() * 2);

	if(!SILENT)
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

	this->literal_lines.push_back(std::string(&stream[0],stream.size()-1));
	stream.clear();
	return true;
}

bool VCFHeader::__parseFirstLine(const char* const data, U32& offset){
	offset = 4;
	if(strncmp(&data[offset], &VCF::Constants::HEADER_VCF_FORMAT[0], VCF::Constants::HEADER_VCF_FORMAT.size()) != 0){
		std::cerr << Helpers::timestamp("ERROR", "BCF") << "Invalid VCF format..." << std::endl;
		std::cerr << std::string(&data[offset], 100) << std::endl;
		this->error_bit = VCF_ERROR_LINE1;
		return false;
	}

	if(strncmp(&data[offset], &VCF::Constants::HEADER_VCF_VERSION[0], VCF::Constants::HEADER_VCF_VERSION.size()) != 0){
		std::cerr << Helpers::timestamp("ERROR", "BCF") << "Invalid VCF version < 4.x..." << std::endl;
		this->error_bit = VCF_ERROR_LINE1;
		return false;
	}

	const char* hit = std::strchr(&data[offset], '\n');
	this->literal_lines.push_back(std::string(&data[offset], hit - &data[offset]));
	offset += (hit + 1) - &data[offset];
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

bool VCFHeader::__parseHeaderLines(const char* const data, U32& offset){
	std::istringstream is(&data[offset]);
	std::string line;
	while(std::getline(is, line)){
		if(line[1] != '#')
			break;

		if(!this->checkLine(&line[0], line.size() - 1)){
			std::cerr << Helpers::timestamp("ERROR", "BCF") << "Failed to validate header lines" << std::endl;
			this->error_bit = VCF_ERROR_LINES;
			return false;
		}
	}

	offset += (U32)is.tellg() - line.size() - 1;

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
		char* found = std::find(&stream[search_position], &stream[stream.size()-1], Tomahawk::VCF::Constants::VCF_DELIMITER);
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
	S32* retValue;
	char* found = 0;
	while(found != &stream[stream.size()-1]){ // while there are samples in line
		found = std::find(&stream[search_position], &stream[stream.size()-1], Tomahawk::VCF::Constants::VCF_DELIMITER);
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

bool VCFHeader::__parseSampleLine(const char* const data, U32& offset, const U32& length){
	// At broken position is main header line
	// Validate header
	if(strncmp(&Tomahawk::VCF::Constants::HEADER_COLUMN[0], &data[offset], Tomahawk::VCF::Constants::HEADER_COLUMN.size()) != 0){
		std::cerr << Helpers::timestamp("ERROR", "VCF") << "Could not validate header line" << std::endl;
		this->error_bit = VCF_ERROR_SAMPLE;
		return false;
	}

	offset += Tomahawk::VCF::Constants::HEADER_COLUMN.size() + 1;
	U64 delimiters_found = 0;
	U32 offset_original = offset;

	while(true){ // while there is samples in line
		const char* const found = std::strchr(&data[offset], Tomahawk::VCF::Constants::VCF_DELIMITER);
		//std::cerr << (void*)found << '\t' << (void*)&data[length] << std::endl;
		if(found == 0 || (*found != Tomahawk::VCF::Constants::VCF_DELIMITER)){
			std::string sampleName(&data[offset], (&data[length - 1] - &data[offset]) - 1); // -2 because offset is +1 and newline is +1
			//std::cerr << sampleName << std::endl;
			++delimiters_found;
			break;
		}

		std::string sampleName(&data[offset], (found - &data[offset]));
		if(sampleName == "FORMAT"){
			offset += found - &data[offset] + 1;
			continue;
		}

		offset += found - &data[offset] + 1;
		++delimiters_found;
	}

	this->buildSampleTable(delimiters_found);

	offset = offset_original;
	S32* retValue;
	while(true){ // while there is samples in line
		const char* const found = std::strchr(&data[offset], Tomahawk::VCF::Constants::VCF_DELIMITER);
		if(found == 0 || (*found != Tomahawk::VCF::Constants::VCF_DELIMITER)){
			std::string sampleName(&data[offset], (&data[length - 1] - &data[offset]) - 1); // -2 because offset is +1 and newline is +1
			if(!this->getSample(sampleName, retValue))
				this->addSample(sampleName);
			else {
				std::cerr << Helpers::timestamp("ERROR", "VCF") << "Duplicated sample name in header..." << std::endl;
				this->error_bit = VCF_ERROR_LINES;
			}
			break;
		}


		std::string sampleName(&data[offset], (found - &data[offset]));
		if(sampleName == "FORMAT"){
			offset += found - &data[offset] + 1;
			continue;
		}

		if(!this->getSample(sampleName, retValue))
			this->addSample(sampleName);
		else {
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Duplicated sample name in header..." << std::endl;
			this->error_bit = VCF_ERROR_LINES;
		}
		offset += found - &data[offset] + 1;
	}

	return true;
}


}
} /* namespace Tomahawk */
