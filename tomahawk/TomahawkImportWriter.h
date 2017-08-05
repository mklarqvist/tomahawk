#ifndef TOMAHAWKIMPORTWRITER_H_
#define TOMAHAWKIMPORTWRITER_H_

#include <fstream>

#include "../TypeDefinitions.h"
#include "../io/vcf/VCFHeaderConstants.h"
#include "../io/vcf/VCFLines.h"
#include "../io/vcf/VCFHeader.h"
#include "../io/BasicBuffer.h"
#include "../io/BasicWriters.h"
#include "../io/GZController.h"
#include "../totempole/TotempoleEntry.h"
#include "TomahawkEntryMeta.h"
#include "../algorithm/compression/TomahawkImportRLE.h"
#include "../algorithm/compression/ByteReshuffle.h"

#include "../totempole/TotempoleReader.h"

namespace Tomahawk {

class TomahawkImportWriterTemp {
	typedef TomahawkImportWriterTemp self_type;
	typedef Tomahawk::IO::WriterFile writer_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::GZController tgzf_type;
	typedef VCF::VCFHeader header_type;
	typedef TotempoleEntry totempole_type;

public:
	void add(buffer_type& rle, buffer_type& meta, totempole_type& totempole_entry){
		totempole_entry.byte_offset = this->totempole_output.getStream().tellp(); // IO offset in Tomahawk output
		this->tgzf_controller.Deflate(meta, rle); // Deflate block
		this->tomahawk_output.getNativeStream() << this->tgzf_controller; // Write tomahawk output
		this->tgzf_controller.Clear();

		// Keep track of largest block observed
		if(meta.size() > this->largest_uncompressed_block)
			this->largest_uncompressed_block = meta.size();

		totempole_entry.uncompressed_size = meta.size(); // Store uncompressed size
		this->totempole_output.getNativeStream() << totempole_entry; // Write totempole output
		++this->blocksWritten; // update number of blocks written
		this->variants_written += totempole_entry.variants; // update number of variants written
	}

	void WriteHeaders(const header_type& header){
		std::ofstream& totempole = this->totempole_output.getNativeStream();
		std::ofstream& tomahawk  = this->tomahawk_output.getNativeStream();

		// Write MAGIC
		this->totempole_output.write(Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH);
		this->tomahawk_output.write(Constants::WRITE_HEADER_MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH);

		// Write samples
		const U64& samples = header.samples;
		Totempole::TotempoleHeader h(samples);
		totempole << h;
		Totempole::TotempoleHeaderBase* hB = reinterpret_cast<Totempole::TotempoleHeaderBase*>(&h);
		tomahawk << *hB;

		// Write out dummy variable for IO offset
		U32 nothing = 0; // Dummy variable
		size_t posOffset = totempole.tellp(); // remember current IO position
		this->totempole_output.write(reinterpret_cast<const char*>(&nothing), sizeof(U32)); // data offset

		// Write the number of contigs
		const U32 n_contigs = header.contigs.size();
		this->totempole_output.write(reinterpret_cast<const char*>(&n_contigs), sizeof(U32));


		// Write contig data to Totempole
		// length | n_char | chars[0 .. n_char - 1]
		for(U32 i = 0; i < header.contigs.size(); ++i){
			Totempole::TotempoleContigBase contig(header.contigs[i].length,
												  header.contigs[i].name.size(),
												  header.contigs[i].name);

			totempole << contig;
		}

		// Write sample names
		// n_char | chars[0..n_char - 1]
		for(U32 i = 0; i < samples; ++i){
			const U32 n_char = header.sampleNames[i].size();
			this->totempole_output.write(reinterpret_cast<const char*>(&n_char), sizeof(U32));
			this->totempole_output.write(reinterpret_cast<const char*>(&header.sampleNames[i][0]), n_char);
		}

		U32 curPos = totempole.tellp(); // remember current IO position
		totempole.seekp(posOffset); // seek to previous position
		totempole.write(reinterpret_cast<const char*>(&curPos), sizeof(U32)); // overwrite data offset
		totempole.seekp(curPos); // seek back to current IO position
	}

private:
	bool writeHeaders(header_type& header);

private:
	U32 blocksWritten;				// number of blocks written
	U32 variants_written;			// number of variants written
	U32 largest_uncompressed_block;// size of largest block in b
	writer_type tomahawk_output;
	writer_type totempole_output;
	tgzf_type tgzf_controller;
};

class TomahawkImportWriter {
public:
	TomahawkImportWriter() :
		blocksWritten_(0),
		variants_written_(0),
		largest_uncompressed_block_(0),
		rleController_(nullptr),
		buffer_rle_(Constants::WRITE_BLOCK_SIZE*2),
		buffer_meta_(Constants::WRITE_BLOCK_SIZE*2),
		//buffer_rle2_(Constants::WRITE_BLOCK_SIZE*2),
		//buffer_meta2_(Constants::WRITE_BLOCK_SIZE*2),
		//buffer_debug1(Constants::WRITE_BLOCK_SIZE*2),
		vcf_header_(nullptr)
	{}

	~TomahawkImportWriter(){
		delete this->rleController_;
		this->buffer_rle_.deleteAll();
		this->buffer_meta_.deleteAll();
	}

	bool Open(const std::string output){
		this->filename = output;
		this->DetermineBasePath();
		this->streamTomahawk.open(this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX, std::ios::binary);
		this->streamTotempole.open(this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX + '.' + Constants::OUTPUT_INDEX_SUFFIX, std::ios::binary);

		// Check streams
		if(!this->streamTomahawk.good()){
			std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX << "!" << std::endl;
			return false;
		}
		if(!this->streamTotempole.good()){
			std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX + '.' + Constants::OUTPUT_INDEX_SUFFIX << "!" << std::endl;
			return false;
		}

		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX << "..." << std::endl;
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX + '.' + Constants::OUTPUT_INDEX_SUFFIX << "..." << std::endl;

		// Write Tomahawk and Totempole headers
		this->WriteHeaders();
		return true;
	}

	void WriteHeaders(void){
		if(this->vcf_header_ == nullptr){
			std::cerr << Helpers::timestamp("ERROR", "INTERNAL") << "Header not set" << std::endl;
			exit(1);
		}

		this->streamTotempole.write(Constants::WRITE_HEADER_INDEX_MAGIC, Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH);
		this->streamTomahawk.write(Constants::WRITE_HEADER_MAGIC, Constants::WRITE_HEADER_MAGIC_LENGTH);

		const U64& samples = this->vcf_header_->samples;
		Totempole::TotempoleHeader h(samples);
		this->streamTotempole << h;
		Totempole::TotempoleHeaderBase* hB = reinterpret_cast<Totempole::TotempoleHeaderBase*>(&h);
		this->streamTomahawk << *hB;

		// Write out dummy variable for IO offset
		U32 nothing = 0; // Dummy variable
		size_t posOffset = this->streamTotempole.tellp(); // remember current IO position
		this->streamTotempole.write(reinterpret_cast<const char*>(&nothing), sizeof(U32)); // data offset

		// Write the number of contigs
		const U32 n_contigs = this->vcf_header_->contigs.size();
		this->streamTotempole.write(reinterpret_cast<const char*>(&n_contigs), sizeof(U32));


		// Write contig data to Totempole
		// length | n_char | chars[0 .. n_char - 1]
		for(U32 i = 0; i < this->vcf_header_->contigs.size(); ++i){
			Totempole::TotempoleContigBase contig(this->vcf_header_->contigs[i].length,
											      this->vcf_header_->contigs[i].name.size(),
											      this->vcf_header_->contigs[i].name);

			this->streamTotempole << contig;
		}

		// Write sample names
		// n_char | chars[0..n_char - 1]
		for(U32 i = 0; i < samples; ++i){
			//const U32 n_char = this->vcf_header_->sampleNames[i].size();
			//this->streamTotempole.write(reinterpret_cast<const char*>(&n_char), sizeof(U32));
			//this->streamTotempole.write(reinterpret_cast<const char*>(&this->vcf_header_->sampleNames[i][0]), n_char);
		}

		U32 curPos = this->streamTotempole.tellp(); // remember current IO position
		this->streamTotempole.seekp(posOffset); // seek to previous position
		this->streamTotempole.write(reinterpret_cast<const char*>(&curPos), sizeof(U32)); // overwrite data offset
		this->streamTotempole.seekp(curPos); // seek back to current IO position
	}

	void WriteFinal(void){
		// Write EOF
		for(U32 i = 0; i < Constants::eof_length; ++i){
			this->streamTotempole.write(reinterpret_cast<const char*>(&Constants::eof[i]), sizeof(U64));
			this->streamTomahawk.write(reinterpret_cast<const char*>(&Constants::eof[i]), sizeof(U64));
		}

		// Re-open file and overwrite block counts and offset
		const U32 shift = Constants::WRITE_HEADER_MAGIC_INDEX_LENGTH + sizeof(float) + sizeof(U64) + sizeof(BYTE);
		this->streamTotempole.flush();
		std::fstream streamTemp(this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX + '.' + Constants::OUTPUT_INDEX_SUFFIX, std::ios_base::binary | std::ios_base::out | std::ios_base::in);

		if(!streamTemp.good()){
			std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not re-open file!" << std::endl;
			exit(1);
		}

		streamTemp.seekg(shift);
		streamTemp.write(reinterpret_cast<const char*>(&this->blocksWritten_), sizeof(U32));
		streamTemp.write(reinterpret_cast<const char*>(&this->largest_uncompressed_block_), sizeof(U32));
		streamTemp.flush();
		streamTemp.close();
	}

	void setHeader(VCF::VCFHeader& header){
		this->vcf_header_ = &header;
		this->rleController_ = new Algorithm::TomahawkImportRLE(header);
		this->rleController_->DetermineBitWidth();
	}

	inline void reset(void){
		this->buffer_rle_.reset();
		this->buffer_meta_.reset();
	}

	inline void operator+=(const VCF::VCFLine& line){
		if(this->totempole_entry_.minPosition == 0)
			this->totempole_entry_.minPosition = line.position;

		this->totempole_entry_.maxPosition = line.position;
		++this->totempole_entry_;
		this->buffer_meta_ += line.position;
		this->buffer_meta_ += line.ref_alt;
		this->rleController_->RunLengthEncode(line, this->buffer_meta_, this->buffer_rle_);
	}

	inline void TotempoleSwitch(const U32 contig, const U32 minPos){
		this->totempole_entry_.reset();
		this->totempole_entry_.contigID = contig;
		this->totempole_entry_.minPosition = minPos;
	}

	// flush and write
	bool flush(void){
		if(this->buffer_meta_.size() == 0){
			std::cerr << Helpers::timestamp("ERROR", "Writer") << "Cannot flush writer with 0 entries..." << std::endl;
			return false;
		}

		this->totempole_entry_.byte_offset = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
		this->gzip_controller_.Deflate(this->buffer_meta_, this->buffer_rle_); // Deflate block
		this->streamTomahawk << this->gzip_controller_; // Write tomahawk output
		this->gzip_controller_.Clear(); // Clean up gzip controller

		// Keep track of largest block observed
		if(this->buffer_meta_.size() > this->largest_uncompressed_block_)
			this->largest_uncompressed_block_ = this->buffer_meta_.size();

		this->totempole_entry_.uncompressed_size = this->buffer_meta_.size(); // Store uncompressed size
		this->streamTotempole << this->totempole_entry_; // Write totempole output
		++this->blocksWritten_; // update number of blocks written
		this->variants_written_ += this->totempole_entry_.variants; // update number of variants written

		this->reset(); // reset buffers
		return true;
	}

	inline bool checkSize() const{
		// if the current size is larger than our desired output block size, return TRUE to trigger a flush
		if(this->buffer_rle_.size() >= Constants::WRITE_BLOCK_SIZE)
			return true;

		return false;
	}

	const U32& blocksWritten(void) const{ return this->blocksWritten_; }
	const U64& size(void) const{ return this->buffer_rle_.size(); }

	void DetermineBaseName(const std::string& string){
		if(string.size() == 0)
			return;

		std::vector<std::string> data = Tomahawk::Helpers::split(string, '.');
		if(data.size() == 1){
			this->baseName = string;
		} else {
			std::transform(data.back().begin(), data.back().end(), data.back().begin(), ::tolower);
			// Data already terminates in the correct output suffix
			if(data.back() == Constants::OUTPUT_SUFFIX){
				this->baseName = string.substr(0, string.find_last_of('.'));
				return;
			}
			this->baseName = string;
		}
	}

	void DetermineBasePath(void){
		if(this->filename.size() == 0)
			return;

		std::vector<std::string> data = Tomahawk::Helpers::split(this->filename, '/');
		if(data.size() == 1){ // If there is no path then we only need to set the basename;
			this->basePath = "";
			this->DetermineBaseName(this->filename);
		} else {
			const U32 lastPos = this->filename.find_last_of('/') + 1;
			this->basePath = this->filename.substr(0, lastPos);
			this->DetermineBaseName(this->filename.substr(lastPos, this->filename.length()));
		}
	}

	U32 GetVariantsWritten(void) const{ return this->variants_written_; }
	TotempoleEntry& getTotempoleEntry(void){ return(this->totempole_entry_); }

private:
	std::ofstream streamTomahawk;	// stream
	std::ofstream streamTotempole;	// stream
	U32 blocksWritten_;				// number of blocks written
	U32 variants_written_;			// number of variants written
	U32 largest_uncompressed_block_;// size of largest block in b

	TotempoleEntry totempole_entry_;
	IO::GZController gzip_controller_;				// todo: writer should have no knowledge of this
	Algorithm::TomahawkImportRLE* rleController_;	// todo: writer should have no knowledge of this
	IO::BasicBuffer buffer_rle_;	// run lengths
	IO::BasicBuffer buffer_meta_;	// meta data for run lengths (chromosome, position, ref/alt)

	VCF::VCFHeader* vcf_header_;

	std::string filename;
	std::string basePath;
	std::string baseName;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWKIMPORTWRITER_H_ */
