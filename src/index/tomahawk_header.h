#ifndef INDEX_TOMAHAWK_HEADER_H_
#define INDEX_TOMAHAWK_HEADER_H_

#include <cassert>

#include "index_contig.h"
#include "../io/BasicBuffer.h"
#include "../algorithm/OpenHashTable.h"
#include "../tomahawk/tomahawk_magic_header.h"
#include "../io/compression/TGZFController.h"

namespace Tomahawk{

/**<
 * This container handles the header data for
 * a `twk` file
 */
class TomahawkHeader{
private:
	typedef TomahawkHeader                    self_type;
	typedef Totempole::HeaderContig           contig_type;
    typedef Base::TomahawkMagicHeader         magic_type;
    typedef IO::BasicBuffer                   buffer_type;
    typedef Hash::HashTable<std::string, S32> hash_table;
    typedef IO::TGZFController                compressor_type;

public:
    TomahawkHeader(void) :
    	contigs_(nullptr),
		sample_names_(nullptr),
		contigs_hash_table_(nullptr),
		sample_hash_table_(nullptr)
	{

	}

    TomahawkHeader(const buffer_type& buffer);
    TomahawkHeader(const char* const data, const U32 l_data);
    TomahawkHeader(const self_type& other);

    // Standard dtor
    ~TomahawkHeader(void){
    	delete [] this->contigs_;
    	delete [] this->sample_names_;
    	delete this->contigs_hash_table_;
    	delete this->sample_hash_table_;
    }

    // Open and close functions
    int open(std::ifstream& stream = std::cin){
    	if(stream.good() == false){
    		std::cerr << Helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
    		return(-1);
    	}

    	stream >> this->magic_;
    	if(this->validate() == false){
    		std::cerr << Helpers::timestamp("ERROR") << "Failed to validate MAGIC header!" << std::endl;
    		return(-2);
    	}

    	if(stream.good() == false){
			std::cerr << Helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
			return(-1);
		}

    	// Parse literal block
    	compressor_type tgzf_controller(this->magic_.l_header_uncompressed + 1024);
    	buffer_type buffer(this->magic_.l_header + 1024);
    	stream.read(buffer.data(), this->magic_.l_header);

    	if(stream.good() == false){
			std::cerr << Helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
			return(-1);
		}

		if(!tgzf_controller.InflateBlock(stream, buffer)){
			std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to get deflate literal TGZF DATA!" << std::endl;
			return(-3);
		}

		// Parse contigs
		// Parse names
		// Construct hash tables
		U32 buffer_position = 0;

		// Parse contigs
		this->contigs_ = new contig_type[this->magic_.getNumberContigs()];
		for(U32 i = 0; i < this->magic_.getNumberContigs(); ++i){
			buffer_position += this->contigs_[i].interpret(&buffer[buffer_position]);
			assert(buffer_position < buffer.size());
		}

		// Parse sample names
		// Encoded as |length in characters|character buffer|
		this->sample_names_ = new std::string[this->magic_.getNumberSamples()];
		for(U32 i = 0; i < this->magic_.getNumberSamples(); ++i){
			const U32 length = *reinterpret_cast<const U32*>(&buffer[buffer_position]);
			buffer_position += sizeof(U32);

			this->sample_names_[i] = std::string(&buffer[buffer_position], length);
			buffer_position += length;
			assert(buffer_position < buffer.size());
		}

		// Remainder is literal data
		const U32 l_literals = buffer.size() - buffer_position;
		this->literals_ = std::string(&buffer[buffer_position], l_literals);

		// Build hash tables for contigs and sample names
		if(this->BuildHashTables() == false){
			std::cerr << Helpers::timestamp("ERROR") << "Cannot build hash tables" << std::endl;
			return(-4);
		}

		// Buffer cleanup
		buffer.deleteAll();

		return(1);
    }

    int write(std::ostream& stream = std::cout){
    	if(stream.good() == false){
			std::cerr << Helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
			return(-1);
		}

    	// Compute uncompressed size
    	const U32 l_uncompressed_size = this->DetermineUncompressedSize();

		buffer_type buffer(l_uncompressed_size + 1024);
		for(U32 i = 0; i < this->magic_.getNumberContigs(); ++i)
			buffer += this->contigs_[i];

		for(U32 i = 0; i < this->magic_.getNumberSamples(); ++i){
			buffer += (U32)this->sample_names_[i].size();
			buffer.Add(this->sample_names_[i].data(), this->sample_names_[i].size());
		}

		buffer.Add(this->literals_.data(), this->literals_.size());
		this->magic_.l_header_uncompressed = buffer.size();
		assert(this->magic_.l_header_uncompressed == l_uncompressed_size);

		compressor_type tgzf_controller(this->magic_.l_header_uncompressed + 1024);
		if(!tgzf_controller.Deflate(buffer)){
			std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed to get deflate literal TGZF DATA!" << std::endl;
			return(-3);
		}

		// Store compressed size
		this->magic_.l_header = tgzf_controller.buffer.size();

		stream << this->magic_;
		if(stream.good() == false){
			std::cerr << Helpers::timestamp("ERROR") << "Stream is bad!" << std::endl;
			return(-1);
		}
		stream.write(tgzf_controller.buffer.data(), tgzf_controller.buffer.size());

		// Cleanup buffer
		buffer.deleteAll();

		return(1);
    }

    // Accessors
    inline std::string& getLiterals(void){ return(this->literals_); }
	inline const std::string& getLiterals(void) const{ return(this->literals_); }
	inline std::string& getSample(const U32 position){ return(this->sample_names_[position]); }
	inline const std::string& getSample(const U32 position) const{ return(this->sample_names_[position]); }

	inline const bool getSample(const std::string& sample_name, const std::string*& return_target) const{
		if(this->sample_hash_table_ == nullptr)
			return false;

		if(sample_name.size() == 0)
			return false;

		if(this->sample_hash_table_->occupied() == 0)
			return false;

		S32* target = nullptr;
		if(this->sample_hash_table_->GetItem(&sample_name[0], &sample_name, target, sample_name.length())){
			return_target = this->sample_names_[*target];
			return true;
		}
		return false;
	}

	inline const bool getContig(const std::string& contig_name, const contig_type*& return_target) const{
		if(this->contigs_hash_table_ == nullptr)
			return false;

		if(contig_name.size() == 0)
			return false;

		if(this->contigs_hash_table_->occupied() == 0)
			return false;

		S32* target = nullptr;
		if(this->contigs_hash_table_->GetItem(&contig_name[0], &contig_name, target, contig_name.length())){
			return_target = this->contig_type[*target];
			return true;
		}
		return false;
	}

	// Updater
	inline void addLiteral(const std::string& string){ this->literals_ += string; }

private:
	inline const bool validate(void) const{ return(this->magic_.validate()); }

	bool BuildHashTables(void){
		// For contigs
		if(this->magic_.getNumberContigs() * 2 < 1024)
			this->contigs_hash_table_ = new hash_table(1024);
		else
			this->contigs_hash_table_ = new hash_table(this->magic_.getNumberContigs() * 2);

		S32* retValue = 0;
		for(U32 i = 0; i < this->magic_.getNumberContigs(); ++i){
			if(this->contigs_hash_table_->GetItem(&this->contigs_[i].name[0], &this->contigs_[i].name, retValue, this->contigs_[i].name.size())){
				std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Duplicated contig! Impossible!" << std::endl;
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
				std::cerr << Helpers::timestamp("ERROR", "TOTEMPOLE") << "Duplicated name! Impossible!" << std::endl;
				return false;
			}
			this->sample_hash_table_->SetItem(&this->sample_names_[i][0], &this->sample_names_[i], i, this->sample_names_[i].size());
		}

		return true;
	}

	const U32 DetermineUncompressedSize(void) const{
		// Encoding of strings
		U32 l_uncompressed = this->magic_.getNumberSamples() * sizeof(U32);
		for(U32 i = 0; i < this->magic_.getNumberSamples(); ++i)
			l_uncompressed += this->sample_names_->size();

		for(U32 i = 0; i < this->magic_.getNumberContigs(); ++i)
			l_uncompressed += this->contigs_[i].name.size() + sizeof(U32) + sizeof(U32);

		l_uncompressed += this->literals_.size();
		return(l_uncompressed);
	}

public:
	magic_type     magic_;    // magic header
	std::string    literals_; // literal data
	contig_type*   contigs_;  // contig data
	std::string*   sample_names_;  // sample names
	hash_table*    contigs_hash_table_; // contig name hash table
	hash_table*    sample_hash_table_;  // sample name hash table
};

}



#endif /* INDEX_TOMAHAWK_HEADER_H_ */
