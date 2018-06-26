#include "tomahawk_header.h"

namespace tomahawk{

TomahawkHeader::TomahawkHeader(void) :
    contigs_(nullptr),
	sample_names_(nullptr),
	contigs_hash_table_(nullptr),
	sample_hash_table_(nullptr)
{

}

// Standard dtor
TomahawkHeader::~TomahawkHeader(void){
	delete [] this->contigs_;
	delete [] this->sample_names_;
	delete this->contigs_hash_table_;
	delete this->sample_hash_table_;
}

// Open and close functions
int TomahawkHeader::open(std::istream& stream){
	if(stream.good() == false){
		std::cerr << helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
		return(-1);
	}

	stream >> this->magic_;
	if(this->validate() == false){
		std::cerr << helpers::timestamp("ERROR") << "Failed to validate MAGIC header!" << std::endl;
		return(-2);
	}

	if(stream.good() == false){
			std::cerr << helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
			return(-1);
	}

	// Parse literal block
	compressor_type tgzf_controller(this->magic_.l_header_uncompressed + 1024);
	buffer_type buffer(this->magic_.l_header + 1024);
	buffer_type buffer_uncompressed(this->magic_.l_header_uncompressed + 1024);
	stream.read(buffer.data(), this->magic_.l_header);
	buffer.n_chars = this->magic_.l_header;

	if(stream.good() == false){
			std::cerr << helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
			return(-1);
		}

		if(!tgzf_controller.Inflate(buffer, buffer_uncompressed)){
			std::cerr << helpers::timestamp("ERROR", "TGZF") << "Failed to get deflate literal TGZF DATA!" << std::endl;
			return(-3);
		}

		// Parse contigs
		// Parse names
		// Construct hash tables
		U32 buffer_position = 0;

		// Parse contigs
		this->contigs_ = new contig_type[this->magic_.getNumberContigs()];
		for(U32 i = 0; i < this->magic_.getNumberContigs(); ++i){
			buffer_position += this->contigs_[i].interpret(&buffer_uncompressed[buffer_position]);
			assert(buffer_position < buffer_uncompressed.size());
		}

		// Parse sample names
		// Encoded as |length in characters|character buffer|
		this->sample_names_ = new std::string[this->magic_.getNumberSamples()];
		for(U32 i = 0; i < this->magic_.getNumberSamples(); ++i){
			const U32 length = *reinterpret_cast<const U32*>(&buffer_uncompressed[buffer_position]);
			buffer_position += sizeof(U32);

			this->sample_names_[i] = std::string(&buffer_uncompressed[buffer_position], length);
			buffer_position += length;
			assert(buffer_position < buffer_uncompressed.size());
		}

		// Remainder is literal data
		const U32 l_literals = buffer_uncompressed.size() - buffer_position;
		this->literals_ = std::string(&buffer_uncompressed[buffer_position], l_literals);

		// Build hash tables for contigs and sample names
		if(this->BuildHashTables() == false){
			std::cerr << helpers::timestamp("ERROR") << "Cannot build hash tables" << std::endl;
			return(-4);
		}

		return(1);
}

int TomahawkHeader::write(std::ostream& stream){
	if(stream.good() == false){
		std::cerr << helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
		return(-1);
	}

	// Compute uncompressed size
	const U32 l_uncompressed_size = this->DetermineUncompressedSize();

	buffer_type buffer(l_uncompressed_size + 1024);
	for(U32 i = 0; i < this->magic_.getNumberContigs(); ++i){
		//std::cerr << helpers::timestamp("DEBUG") << this->contigs_[i] << std::endl;
		buffer += this->contigs_[i];
	}

	for(U32 i = 0; i < this->magic_.getNumberSamples(); ++i){
		buffer += (U32)this->sample_names_[i].size();
		//std::cerr << helpers::timestamp("DEBUG") << this->sample_names_[i] << std::endl;
		buffer.Add(this->sample_names_[i].data(), this->sample_names_[i].size());
	}

	buffer.Add(this->literals_.data(), this->literals_.size());
	this->magic_.l_header_uncompressed = buffer.size();
	//std::cerr << buffer.size() << "\t" << l_uncompressed_size << std::endl;
	assert(this->magic_.l_header_uncompressed == l_uncompressed_size);


	compressor_type tgzf_controller(this->magic_.l_header_uncompressed + 1024);
	if(!tgzf_controller.Deflate(buffer)){
		std::cerr << helpers::timestamp("ERROR", "TGZF") << "Failed to get deflate literal TGZF DATA!" << std::endl;
		return(-3);
	}

	// Store compressed size
	this->magic_.l_header = tgzf_controller.buffer.size();

	stream << this->magic_;
	if(stream.good() == false){
		std::cerr << helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
		return(-1);
	}
	stream.write(tgzf_controller.buffer.data(), tgzf_controller.buffer.size());

	//std::cerr << helpers::timestamp("DEBUG") << this->magic_.l_header << "->" << this->magic_.l_header_uncompressed << '\t' << buffer.size() << "/" << buffer.capacity() << std::endl;

	return(tgzf_controller.buffer.size());
}

const bool TomahawkHeader::getSample(const std::string& sample_name, const std::string*& return_target) const{
	if(this->sample_hash_table_ == nullptr)
		return false;

	if(sample_name.size() == 0)
		return false;

	if(this->sample_hash_table_->occupied() == 0)
		return false;

	S32* target = nullptr;
	if(this->sample_hash_table_->GetItem(&sample_name[0], &sample_name, target, sample_name.length())){
		return_target = &this->sample_names_[*target];
		return true;
	}
	return false;
}

const bool TomahawkHeader::getContigName(const std::string& contig_name, const std::string*& return_target) const{
	if(this->contigs_hash_table_ == nullptr)
		return false;

	if(contig_name.size() == 0)
		return false;

	if(this->contigs_hash_table_->occupied() == 0)
		return false;

	S32* target = nullptr;
	if(this->contigs_hash_table_->GetItem(&contig_name[0], &contig_name, target, contig_name.length())){
		return_target = &this->contigs_[*target].name;
		return true;
	}
	return false;
}

const bool TomahawkHeader::getContig(const std::string& contig_name, const contig_type*& return_target) const{
	if(this->contigs_hash_table_ == nullptr)
		return false;

	if(contig_name.size() == 0)
		return false;

	if(this->contigs_hash_table_->occupied() == 0)
		return false;

	S32* target = nullptr;
	if(this->contigs_hash_table_->GetItem(&contig_name[0], &contig_name, target, contig_name.length())){
		return_target = &this->contigs_[*target];
		return true;
	}
	return false;
}

const S32 TomahawkHeader::getContigID(const std::string& contig_name) const{
	if(this->contigs_hash_table_ == nullptr)
		return false;

	if(contig_name.size() == 0)
		return false;

	if(this->contigs_hash_table_->occupied() == 0)
		return false;

	S32* target = nullptr;
	if(this->contigs_hash_table_->GetItem(&contig_name[0], &contig_name, target, contig_name.length())){
		return(*target);
	}
	return(-1);
}

bool TomahawkHeader::BuildHashTables(void){
	// For contigs
	if(this->magic_.getNumberContigs() * 2 < 1024)
		this->contigs_hash_table_ = new hash_table(1024);
	else
		this->contigs_hash_table_ = new hash_table(this->magic_.getNumberContigs() * 2);

	S32* retValue = 0;
	for(U32 i = 0; i < this->magic_.getNumberContigs(); ++i){
		if(this->contigs_hash_table_->GetItem(&this->contigs_[i].name[0], &this->contigs_[i].name, retValue, this->contigs_[i].name.size())){
			std::cerr << helpers::timestamp("ERROR", "TOTEMPOLE") << "Duplicated contig! Impossible!" << std::endl;
			return false;
		}
		this->contigs_hash_table_->SetItem(&this->contigs_[i].name[0], &this->contigs_[i].name, i, this->contigs_[i].name.size());
	}

	// For sample names
	if(this->magic_.getNumberSamples() * 2 < 1024)
		this->sample_hash_table_ = new hash_table(1024);
	else
		this->sample_hash_table_ = new hash_table(this->magic_.getNumberSamples() * 2);

	retValue = 0;
	for(U32 i = 0; i < this->magic_.getNumberSamples(); ++i){
		if(this->sample_hash_table_->GetItem(&this->sample_names_[i][0], &this->sample_names_[i], retValue, this->sample_names_[i].size())){
			std::cerr << helpers::timestamp("ERROR", "TOTEMPOLE") << "Duplicated name! Impossible!" << std::endl;
			return false;
		}
		this->sample_hash_table_->SetItem(&this->sample_names_[i][0], &this->sample_names_[i], i, this->sample_names_[i].size());
	}

	return true;
}

const U32 TomahawkHeader::DetermineUncompressedSize(void) const{
	U32 l_uncompressed = 0;
	for(U32 i = 0; i < this->magic_.getNumberContigs(); ++i){
		l_uncompressed += this->contigs_[i].name.size() + sizeof(U32)*2;
	}

	for(U32 i = 0; i < this->magic_.getNumberSamples(); ++i){
		l_uncompressed += this->sample_names_[i].size() + sizeof(U32);
	}

	l_uncompressed += this->literals_.size();

	return(l_uncompressed);
}


}
