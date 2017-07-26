#ifndef TOMAHAWKOUTPUTREADER_H_
#define TOMAHAWKOUTPUTREADER_H_

#include <vector>
#include <fstream>
#include <string>
#include <stddef.h>
#include <regex>

#include "../../TypeDefinitions.h"
#include "../../tomahawk/MagicConstants.h"
#include "TomahawkOutputEntry.h"
#include "../../io/PackedEntryReader.h"
#include "TomahawkOutputFilterController.h"
#include "../../io/BasicBuffer.h"
#include "../../io/GZController.h"
#include "../../io/BasicWriters.h"
#include "../../io/totempole/TotempoleMagic.h"
#include "../../io/TGZFHeader.h"
#include "../../third_party/intervalTree.h"


namespace Tomahawk {
namespace IO {
// Todo: TomahawkOutputIndexReader

class TomahawkOutputReader {
	typedef TomahawkOutputEntry entry_type;
	typedef TomahawkOutputFilterController filter_type;
	typedef PackedEntryReader<entry_type, sizeof(entry_type)> reader_type;
	typedef IO::GenericWriterInterace writer_type;
	typedef TomahawkOutputHeader<Tomahawk::Constants::WRITE_HEADER_LD_MAGIC_LENGTH> header_type;
	typedef Totempole::TotempoleContigBase contig_type;
	typedef TGZFHeader tgzf_type;
	typedef Hash::HashTable<std::string, U32> hash_table;
	typedef IO::GZController tgzf_controller_type;
	typedef Tomahawk::Algorithm::ContigInterval interval_type;
	typedef Tomahawk::Algorithm::IntervalTree<interval_type, U32> tree_type;

public:
	TomahawkOutputReader();
	~TomahawkOutputReader(){
		delete [] this->contigs;
		delete contig_htable;
		if(interval_tree != nullptr){
			for(U32 i = 0; i < this->header.n_contig; ++i)
				delete this->interval_tree[i];
		}
		delete [] interval_tree_entries;

		delete interval_tree;
		delete writer;

		this->buffer.deleteAll();
		this->output_buffer.deleteAll();
	}

	// Streaming functions
	bool getBlock(const U32 blockID);
	bool getBlock(std::vector< std::pair<U32, U32> >& pairs);
	bool getBlocks(void);

	bool AddRegions(std::vector<std::string>& positions){
		if(positions.size() == 0)
			return false;

		if(this->interval_tree == nullptr){
			this->interval_tree = new tree_type*[this->header.n_contig];
			for(U32 i = 0; i < this->header.n_contig; ++i)
				this->interval_tree[i] = nullptr;
		}
		if(this->interval_tree_entries == nullptr)
			this->interval_tree_entries = new std::vector<interval_type>[this->header.n_contig];

		for(U32 i = 0; i < positions.size(); ++i){
			std::cerr << i << ": " << positions[i] << std::endl;
			// Pattern cA:pAf-pAt;cB:pBf-pBt
			if(positions[i].find(',') != std::string::npos){
				std::cerr << "linked intervals" << std::endl;
				std::vector<std::string> ret = Helpers::split(positions[i], ',');
				if(ret.size() == 1){
					std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
					return false;

				} else if(ret.size() == 2){
					// parse left
					interval_type intervalLeft;
					if(!__ParseRegion(ret[0], intervalLeft))
						return false;

					// parse right
					interval_type intervalRight;
					if(!__ParseRegion(ret[1], intervalRight))
						return false;

					// Todo: WARNING
					// This results in illegal pointers if the vector resizes
					// and pointers change
					this->interval_tree_entries[intervalLeft.contigID].push_back(interval_type(intervalLeft));
					this->interval_tree_entries[intervalRight.contigID].push_back(interval_type(intervalRight));
					if(intervalLeft.contigID != intervalRight.contigID){
						this->interval_tree_entries[intervalLeft.contigID].back().value = &this->interval_tree_entries[intervalRight.contigID].back();
						this->interval_tree_entries[intervalRight.contigID].back().value = &this->interval_tree_entries[intervalLeft.contigID].back();

						// Link the intervals together
						std::cerr << this->interval_tree_entries[intervalLeft.contigID].back() << '\t' << this->interval_tree_entries[intervalRight.contigID].back() << std::endl;
					} else {
						this->interval_tree_entries[intervalLeft.contigID].back().value = &this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2];
						this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2].value = &this->interval_tree_entries[intervalLeft.contigID].back();

						std::cerr << this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2] << '\t' << this->interval_tree_entries[intervalRight.contigID].back() << std::endl;
					}

				} else {
					std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << positions[i] << "!" << std::endl;
					return false;
				}
			} else {
				interval_type interval;
				if(!__ParseRegion(positions[i], interval))
					return false;

				this->interval_tree_entries[interval.contigID].push_back(interval_type(interval));
			}
		}

		for(U32 i = 0; i < this->header.n_contig; ++i){
			if(this->interval_tree_entries[i].size() != 0){
				std::cerr << "constructing itree for id: " << i << std::endl;
				this->interval_tree[i] = new tree_type(this->interval_tree_entries[i]);
			} else
				this->interval_tree[i] = nullptr;
		}

		return true;
	}

	bool __ParseRegion(const std::string& region, interval_type& interval){
		std::vector<std::string> ret = Helpers::split(region, ':');
		if(ret.size() == 1){
			if(ret[0].find('-') != std::string::npos){
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
				return false;
			}

			// is contigID only
			std::cerr << "contigONly" << std::endl;
			U32* contigID;
			if(!this->contig_htable->GetItem(&region[0], &region, contigID, region.size())){
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << region << " is not defined in the header!" << std::endl;
				return false;
			}
			interval(0, this->contigs[*contigID].bases, *contigID);

		} else if(ret.size() == 2){
			// is contigID:pos-pos
			U32* contigID;
			if(!this->contig_htable->GetItem(&ret[0][0], &ret[0], contigID, ret[0].size())){
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Contig: " << ret[0] << " is not defined in the header!" << std::endl;
				return false;
			}

			std::vector<std::string> retPos = Helpers::split(ret[1], '-');
			if(retPos.size() == 1){
				// only one pos
				const double pos = std::stod(retPos[0]);
				std::cerr << "single position: " << pos << std::endl;
				interval(pos, pos, *contigID);

			} else if(retPos.size() == 2){
				// is two positions
				double posA = std::stod(retPos[0]);
				double posB = std::stod(retPos[1]);
				if(posB < posA){
					std::cerr << "end position > start position: swapping" << std::endl;
					std::swap(posA, posB);
				}
				std::cerr << "full region: " << this->contigs[*contigID].name << ":" << posA << '-' << posB << std::endl;
				interval(posA, posB, *contigID);

			} else {
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
				return false;
			}
		} else {
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}

		return true;
	}

	bool Open(const std::string input){
		this->stream.open(input, std::ios::binary | std::ios::in | std::ios::ate);
		if(!this->stream.good()){
			std::cerr << "failed to open file " << input << std::endl;
			return false;
		}

		this->filesize = this->stream.tellg();
		this->stream.seekg(0);

		if(!this->stream.good()){
			std::cerr << "bad stream" << std::endl;
			return false;
		}

		this->stream >> this->header;
		if(!this->header.validate(Tomahawk::Constants::WRITE_HEADER_LD_MAGIC)){
			std::cerr << "failed to validate header" << std::endl;
			return false;
		}

		if(!this->ParseHeader()){
			std::cerr << "failed to parse header" << std::endl;
			return false;
		}

		return true;
	}

	bool ParseHeader(void){
		if(this->header.n_contig == 0)
			return false;

		if(this->header.n_contig < 1024)
			this->contig_htable = new hash_table(1024);
		else
			this->contig_htable = new hash_table(this->header.n_contig * 2);

		this->contigs = new contig_type[this->header.n_contig];
		U32* ret;
		for(U32 i = 0; i < this->header.n_contig; ++i){
			this->stream >> this->contigs[i];
			// std::cerr << this->contigs[i] << std::endl;
			if(!this->contig_htable->GetItem(&this->contigs[i].name[0], &this->contigs[i].name, ret, this->contigs[i].name.size())){
				// Add to hash table
				this->contig_htable->SetItem(&this->contigs[i].name[0], &this->contigs[i].name, i, this->contigs[i].name.size());
			} else {
				std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Duplicated contig name: " << this->contigs[i].name << "!" << std::endl;
				exit(1); // unrecoverable error
			}
		}

		return true;
	}

	bool nextBlock(void){
		// Stream died
		if(!this->stream.good()){
			std::cerr << "stream died" << std::endl;
			return false;
		}

		// EOF
		if(this->stream.tellg() == this->filesize){
			std::cerr << "eof" << std::endl;
			return false;
		}

		buffer.resize(sizeof(tgzf_type));
		this->stream.read(&buffer.data[0],  Constants::TGZF_BLOCK_HEADER_LENGTH);
		const tgzf_type* h = reinterpret_cast<const tgzf_type*>(&buffer.data[0]);
		buffer.pointer = Constants::TGZF_BLOCK_HEADER_LENGTH;
		if(!h->Validate()){
			std::cerr << "failed to validate" << std::endl;
			return false;
		}

		buffer.resize(h->BSIZE);

		// Recast because if buffer is resized then the pointer address is incorrect
		// resulting in segfault
		h = reinterpret_cast<const tgzf_type*>(&buffer.data[0]);

		this->stream.read(&buffer.data[Constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - Constants::TGZF_BLOCK_HEADER_LENGTH);
		if(!this->stream.good()){
			std::cerr << "truncated file" << std::endl;
			return false;
		}

		buffer.pointer = h->BSIZE;
		const U32* outsize = reinterpret_cast<const U32*>(&buffer[buffer.pointer -  sizeof(U32)]);
		//const U32* crc = reinterpret_cast<const U32*>(&buffer[buffer.pointer -  sizeof(U32) - sizeof(U32)]);
		//std::cerr << *outsize << '\t' << *crc << std::endl;
		output_buffer.resize(*outsize);
		this->output_buffer.reset();

		if(!this->gzip_controller.Inflate(buffer, output_buffer)){
			std::cerr << "failed inflate" << std::endl;
			return false;
		}

		if(this->output_buffer.size() == 0){
			std::cerr << "empty data" << std::endl;
			return false;
		}

		// Reset buffer
		this->buffer.reset();

		// Reset iterator position and size
		this->position = 0;
		this->size = this->output_buffer.size() / sizeof(TomahawkOutputEntry);

		// Validity check
		if(this->output_buffer.size() % sizeof(TomahawkOutputEntry) != 0){
			std::cerr << "data is corrupted" << std::endl;
			return false;
		}

		return true;
	}

	const TomahawkOutputEntry* operator[](const U32 p) const{ return(reinterpret_cast<TomahawkOutputEntry*>(&this->output_buffer.data[sizeof(TomahawkOutputEntry)*p])); }

	bool nextVariant(const TomahawkOutputEntry*& entry){
		if(this->position == this->size){
			if(!this->nextBlock())
				return false;
		}

		entry = (*this)[this->position];
		++this->position;

		return true;
	}

	// Other
	bool view(const std::string& filename);
	bool index(const std::string& filename);
	bool summary(const std::string& input);

	// Read entire file into memory
	filter_type& getFilter(void){ return this->filter; }

private:
	bool __viewOnly(void);
	bool __viewFilter(void);
	bool __viewRegion(void);

public:
	//U64 samples; 	// has to match header
	//float version;	// has to match header
	U64 filesize;	// input file size

	U32 position;
	U32 size;

	std::ifstream stream; // reader stream
	header_type header; // header

	IO::BasicBuffer buffer; // internal buffer
	IO::BasicBuffer output_buffer; // internal buffer
	tgzf_controller_type gzip_controller; // TGZF controller
	filter_type filter;	// filter parameters
	writer_type* writer; // writer interface
	// Todo: PackedEntryIterator taking as input char* and length or IO::BasicBuffer
	contig_type* contigs;
	hash_table* contig_htable; // map input string to internal contigID
	tree_type** interval_tree;
	std::vector<interval_type>* interval_tree_entries;

	//temp
	reader_type reader;
};

}
} /* namespace Tomahawk */

#endif /* TOMAHAWKOUTPUTREADER_H_ */
