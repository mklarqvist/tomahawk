#include <algorithm>
#include <bitset>
#include <queue>

#include "../../support/helpers.h"
#include "../../support/MagicConstants.h"
#include "../../algorithm/OpenHashTable.h"
#include "TomahawkOutputReader.h"
#include "../../algorithm/sort/TomahawkOutputSort.h"
#include "TomahawkOutputWriter.h"
#include "TomahawkOutputStats.h"

namespace Tomahawk {
namespace IO {

TomahawkOutputReader::TomahawkOutputReader() :
		filesize(0),
		position(0),
		size(0),
		hasIndex(false),
		output_header(true),
		writer_output_type(WRITER_TYPE::natural),
		writer(nullptr),
		contigs(nullptr),
		contig_htable(nullptr),
		interval_tree(nullptr),
		interval_tree_entries(nullptr)
{}

TomahawkOutputReader::~TomahawkOutputReader(){
	delete [] this->contigs;
	delete contig_htable;
	if(interval_tree != nullptr){
		for(U32 i = 0; i < this->header.n_contig; ++i)
			delete this->interval_tree[i];
	}
	//delete interval_tree;
	delete [] interval_tree_entries;
	delete interval_tree;
	this->buffer.deleteAll();
	this->output_buffer.deleteAll();
	delete this->writer;
}

bool TomahawkOutputReader::view(const std::string& input){
	if(this->interval_tree != nullptr) // If regions have been set: use region-filter function
		return(this->__viewRegion());
	else if(this->filter.any_filter_user_set){
		return(this->__viewFilter()); // Otherwise normal filter function
	} else
		return(this->__viewOnly());
}

bool TomahawkOutputReader::OpenWriter(void){
	if(this->writer_output_type == WRITER_TYPE::natural){
		this->writer = new TomahawkOutputWriterNatural(this->contigs, &this->header);
	}
	else this->writer = new TomahawkOutputWriter(this->contigs, &this->header);

	if(!this->writer->open())
		return false;

	if(this->output_header)
		this->writer->writeHeader(this->literals);

	return true;
}

bool TomahawkOutputReader::OpenWriter(const std::string output_file){
	if(this->writer_output_type == WRITER_TYPE::natural){
		this->writer = new TomahawkOutputWriterNatural(this->contigs, &this->header);
	}
	else this->writer = new TomahawkOutputWriter(this->contigs, &this->header);

	if(!this->writer->open(output_file)){
		std::cerr << Helpers::timestamp("ERROR","WRITER") << "Failed to open output file: " << output_file << std::endl;
		return false;
	}

	if(this->output_header)
		this->writer->writeHeader(this->literals);

	return true;
}

bool TomahawkOutputReader::__viewRegion(void){
	if(!this->OpenWriter())
		return false;

	if(this->toi_reader.ERROR_STATE == toi_reader_type::TOI_OK){
		if(this->toi_reader.getIsSortedExpanded()){
			std::cerr << "has sorted and expanded index" << std::endl;
		} else
			std::cerr << "has partial index only" << std::endl;
	} else {
		std::cerr << "no index" << std::endl;
	}

	if(this->interval_tree != nullptr){
		const entry_type*  entry;

		while(this->nextVariant(entry)){
			this->__checkRegionNoIndex(entry);
		} // end while next variant
	}

	return true;
}

bool TomahawkOutputReader::__checkRegionNoIndex(const entry_type* const entry){
	// If iTree for contigA exists
	if(this->interval_tree[entry->AcontigID] != nullptr){
		std::vector<interval_type> rets = this->interval_tree[entry->AcontigID]->findOverlapping(entry->Aposition, entry->Aposition);
		if(rets.size() > 0){
			for(U32 i = 0; i < rets.size(); ++i){
				if(rets[i].value != nullptr){ // if linked
					if((entry->BcontigID == rets[i].value->contigID) &&
					   (entry->Bposition >= rets[i].value->start && entry->Bposition <= rets[i].value->stop)){
						if(this->filter.filter(*entry))
							//entry->write(std::cout, this->contigs);
							*this->writer << entry;

						return true;
					} // end match
				} else { //  not linked
					if(this->filter.filter(*entry))
						//entry->write(std::cout, this->contigs);
						*this->writer << entry;

					return true;
				}
			}
		}
	}

	// If iTree for contigB exists
	if(this->interval_tree[entry->BcontigID] != nullptr){
		std::vector<interval_type> rets = this->interval_tree[entry->BcontigID]->findOverlapping(entry->Bposition, entry->Bposition);
		if(rets.size() > 0){
			for(U32 i = 0; i < rets.size(); ++i){
				if(rets[i].value != nullptr){ // if linked
					if((entry->AcontigID == rets[i].value->contigID) &&
					   (entry->Aposition >= rets[i].value->start && entry->Aposition <= rets[i].value->stop)){
						if(this->filter.filter(*entry)){
							//entry->write(std::cout, this->contigs);
							*this->writer << entry;
						}
						return true;
					} // end match
				} else { // not linked
					if(this->filter.filter(*entry))
						//entry->write(std::cout, this->contigs);
						*this->writer << entry;

					return true;
				}
			}
		} // end if any hit in iTree b
	} // end iTree b

	return false;
}

bool TomahawkOutputReader::__viewOnly(void){
	this->literals += "\n##tomahawk_viewCommand=" + Helpers::program_string();
	this->literals += "\n##tomahawk_viewFilters=" + this->filter.getInterpretedString() + " filter=NO regions=FALSE";

	if(!this->OpenWriter())
		return false;

	// Natural output required parsing
	if(this->writer_output_type == WRITER_TYPE::natural){
		const entry_type* entry;
		while(this->nextVariant(entry))
			*this->writer << entry;

	}
	// Binary output without filtering simply writes it back out
	else if(this->writer_output_type == WRITER_TYPE::binary){
		while(this->nextBlock())
			this->writer->write(this->output_buffer);
	}

	return true;
}

bool TomahawkOutputReader::__viewFilter(void){
	this->literals += "\n##tomahawk_viewCommand=" + Helpers::program_string();
	this->literals += "\n##tomahawk_viewFilters=" + this->filter.getInterpretedString() + " filter=YES regions=FALSE";

	if(!this->OpenWriter())
		return false;

	const entry_type* entry;
	while(this->nextVariant(entry)){
		if(this->filter.filter(*entry))
			*this->writer << entry;
	}

	return true;
}

bool TomahawkOutputReader::AddRegions(std::vector<std::string>& positions){
	if(positions.size() == 0)
		return true;

	if(this->interval_tree == nullptr){
		this->interval_tree = new tree_type*[this->header.n_contig];
		for(U32 i = 0; i < this->header.n_contig; ++i)
			this->interval_tree[i] = nullptr;
	}

	if(this->interval_tree_entries == nullptr)
		this->interval_tree_entries = new std::vector<interval_type>[this->header.n_contig];

	for(U32 i = 0; i < positions.size(); ++i){
		//std::cerr << i << ": " << positions[i] << std::endl;
		// Pattern cA:pAf-pAt;cB:pBf-pBt
		if(positions[i].find(',') != std::string::npos){
			//std::cerr << "linked intervals" << std::endl;
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
					this->interval_tree_entries[intervalLeft.contigID].back().value  = &this->interval_tree_entries[intervalRight.contigID].back();
					this->interval_tree_entries[intervalRight.contigID].back().value = &this->interval_tree_entries[intervalLeft.contigID].back();
				} else {
					this->interval_tree_entries[intervalLeft.contigID].back().value = &this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2];
					this->interval_tree_entries[intervalLeft.contigID][this->interval_tree_entries[intervalLeft.contigID].size() - 2].value = &this->interval_tree_entries[intervalLeft.contigID].back();
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
			if(this->toi_reader.size() > 0){
				this->toi_reader.findOverlap(interval.contigID);
				std::cerr << "has toi index" << std::endl;
				//this->toi_reader.
			}

		}
	}

	for(U32 i = 0; i < this->header.n_contig; ++i){
		if(this->interval_tree_entries[i].size() != 0){
			this->interval_tree[i] = new tree_type(this->interval_tree_entries[i]);
		} else
			this->interval_tree[i] = nullptr;
	}

	return true;
}

bool TomahawkOutputReader::__ParseRegion(const std::string& region, interval_type& interval){
	std::vector<std::string> ret = Helpers::split(region, ':');
	if(ret.size() == 1){
		if(ret[0].find('-') != std::string::npos){
			std::cerr << Helpers::timestamp("ERROR", "INTERVAL") << "Illegal interval: " << region << "!" << std::endl;
			return false;
		}

		// is contigID only
		//std::cerr << "contigONly" << std::endl;
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
			//std::cerr << "single position: " << pos << std::endl;
			interval(pos, pos, *contigID);

		} else if(retPos.size() == 2){
			// is two positions
			double posA = std::stod(retPos[0]);
			double posB = std::stod(retPos[1]);
			if(posB < posA){
				//std::cerr << "end position > start position: swapping" << std::endl;
				std::swap(posA, posB);
			}
			//std::cerr << (U64)posA << '\t' << (U64)posB << std::endl;
			//std::cerr << "full region: " << this->contigs[*contigID].name << ":" << posA << '-' << posB << std::endl;
			interval(posA, posB, *contigID);
			//std::cerr << interval << std::endl;

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

bool TomahawkOutputReader::__Open(const std::string input){
	this->stream.open(input, std::ios::binary | std::ios::in | std::ios::ate);
	if(!this->stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to open file: " << input << std::endl;
		return false;
	}

	if(this->toi_reader.Open(input + "." + Tomahawk::Constants::OUTPUT_LD_SORT_INDEX_SUFFIX)){
		//std::cerr << "could not open index" << std::endl;
		this->hasIndex = true;
	}

	this->filesize = this->stream.tellg();
	this->stream.seekg(0);

	if(!this->stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Bad stream!" << std::endl;
		return false;
	}

	this->stream >> this->header;
	if(!this->header.validate(Tomahawk::Constants::WRITE_HEADER_LD_MAGIC)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to validate header!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::Open(const std::string input){
	if(!this->__Open(input))
		return false;

	if(!this->ParseHeader()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to parse header!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::OpenExtend(const std::string input){
	if(!this->__Open(input))
		return false;

	if(!this->ParseHeaderExtend()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to extend header!" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::__concat(const std::vector<std::string>& files, const std::string& output){
	if(files.size() == 0){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "No input files..." << std::endl;
		return false;
	}

	// open first one
	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "CONCAT") << "Opening input: " << files[0] << "..." << std::endl;

	if(!this->Open(files[0])){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Failed to parse: " << files[0] << "..." << std::endl;
		return false;
	}

	this->setWriterType(0);
	this->setWriteHeader(true);
	this->literals += "\n##tomahawk_concatCommand=" + Helpers::program_string();
	this->literals += "\n##tomahawk_concatFiles=";
	for(U32 i = 0; i < files.size(); ++i)
		this->literals += files[i] + ',';

	if(!this->OpenWriter(output)){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Failed to open writer..." << std::endl;
		return false;
	}

	while(this->nextBlock()){
		this->writer->write(this->output_buffer);
	}

	for(U32 i = 1; i < files.size(); ++i){
		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "CONCAT") << "Opening input: " << files[i] << "..." << std::endl;

		this->stream.close();
		if(!this->OpenExtend(files[i])){
			std::cerr << Helpers::timestamp("ERROR","TWO") << "Failed to parse: " << files[i] << "..." << std::endl;
			return false;
		}

		while(this->nextBlock()){
			this->writer->write(this->output_buffer);
		}
	}

	this->writer->flush();
	this->writer->close();
	return true;
}

bool TomahawkOutputReader::concat(const std::vector<std::string>& files, const std::string& output){
	if(files.size() == 0){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "No input files given..." << std::endl;
		return false;
	}

	return(this->__concat(files, output));
}

bool TomahawkOutputReader::concat(const std::string& file_list, const std::string& output){
	if(file_list.size() == 0){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "No input file list given..." << std::endl;
		return false;
	}

	std::ifstream file_list_read(file_list);
	if(!file_list_read.good()){
		std::cerr << Helpers::timestamp("ERROR","TWO") << "Failed to get file_list..." << std::endl;
		return false;
	}

	std::vector<std::string> files;
	std::string line;
	while(getline(file_list_read, line)){
		if(line.size() == 0){
			std::cerr << Helpers::timestamp("WARNING","TWO") << "Empty line" << std::endl;
			break;
		}
		files.push_back(line);
	}

	return(this->__concat(files, output));
}

bool TomahawkOutputReader::ParseHeader(void){
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

	if(!this->gzip_controller.InflateBlock(this->stream, this->buffer)){
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "Failed to get TWO block" << std::endl;
		return false;
	}

	this->literals = std::string(this->gzip_controller.buffer.data);

	return true;
}

bool TomahawkOutputReader::ParseHeaderExtend(void){
	if(this->header.n_contig == 0)
		return false;

	U32* ret;
	for(U32 i = 0; i < this->header.n_contig; ++i){
		this->stream >> this->contigs[i];
		// std::cerr << this->contigs[i] << std::endl;
		if(!this->contig_htable->GetItem(&this->contigs[i].name[0], &this->contigs[i].name, ret, this->contigs[i].name.size())){
			std::cerr << Helpers::timestamp("ERROR","TWO") << "Contig does not exist in other file" << std::endl;
			return false;
		}
	}

	if(!this->gzip_controller.InflateBlock(this->stream, this->buffer)){
		std::cerr << Helpers::timestamp("ERROR","TGZF") << "Failed to get TWO block" << std::endl;
		return false;
	}

	return true;
}

bool TomahawkOutputReader::nextBlock(void){
	// Stream died
	if(!this->stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Stream died!" << std::endl;
		return false;
	}

	// EOF
	if(this->stream.tellg() == this->filesize){
		//std::cerr << "eof" << std::endl;
		return false;
	}

	buffer.resize(sizeof(tgzf_type));
	this->stream.read(&buffer.data[0],  Constants::TGZF_BLOCK_HEADER_LENGTH);
	const tgzf_type* h = reinterpret_cast<const tgzf_type*>(&buffer.data[0]);
	buffer.pointer = Constants::TGZF_BLOCK_HEADER_LENGTH;
	if(!h->Validate()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to validate!" << std::endl;
		return false;
	}

	buffer.resize(h->BSIZE); // make sure all data will fit

	// Recast because if buffer is (actually) resized then the pointer address is incorrect
	// resulting in segfault
	h = reinterpret_cast<const tgzf_type*>(&buffer.data[0]);

	this->stream.read(&buffer.data[Constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - Constants::TGZF_BLOCK_HEADER_LENGTH);
	if(!this->stream.good()){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Truncated file..." << std::endl;
		return false;
	}

	buffer.pointer = h->BSIZE;
	const U32 uncompressed_size = *reinterpret_cast<const U32*>(&buffer[buffer.pointer -  sizeof(U32)]);
	output_buffer.resize(uncompressed_size);
	this->output_buffer.reset();

	if(!this->gzip_controller.Inflate(buffer, output_buffer)){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed inflate!" << std::endl;
		return false;
	}

	if(this->output_buffer.size() == 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Empty data!" << std::endl;
		return false;
	}

	// Reset buffer
	this->buffer.reset();

	// Reset iterator position and size
	this->position = 0;
	this->size = this->output_buffer.size() / sizeof(entry_type);

	// Validity check
	if(this->output_buffer.size() % sizeof(entry_type) != 0){
		std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Data is corrupted!" << std::endl;
		return false;
	}

	//std::cerr << this->position << '/' << this->size << '\t' << this->output_buffer.size() << '/' << this->output_buffer.capacity() << std::endl;

	return true;
}

bool TomahawkOutputReader::nextBlockUntil(const U32 limit){
	// Check if resize required
	if(this->output_buffer.capacity() < limit + 500536)
		this->output_buffer.resize(limit + 500536);


	this->position = 0;
	this->output_buffer.reset();

	// Keep inflating DATA until bounds is reached
	while(this->output_buffer.size() <= limit){
		//std::cerr << this->stream.tellg() << '/' << this->filesize << std::endl;

		// Stream died
		if(!this->stream.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Stream died!" << std::endl;
			return false;
		}

		// EOF
		// Casting stream to U64 is safe as this point is not
		// reached if above good() return fails
		if((U64)this->stream.tellg() == this->filesize){
			//std::cerr << "eof" << std::endl;
			return false;
		}

		buffer.resize(sizeof(tgzf_type));
		this->stream.read(&buffer.data[0],  Constants::TGZF_BLOCK_HEADER_LENGTH);
		const tgzf_type* h = reinterpret_cast<const tgzf_type*>(&buffer.data[0]);
		buffer.pointer = Constants::TGZF_BLOCK_HEADER_LENGTH;
		if(!h->Validate()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed to validate header!" << std::endl;
			return false;
		}

		buffer.resize(h->BSIZE); // make sure all data will fit

		// Recast because if buffer is resized then the pointer address is incorrect
		// resulting in segfault
		h = reinterpret_cast<const tgzf_type*>(&buffer.data[0]);

		this->stream.read(&buffer.data[Constants::TGZF_BLOCK_HEADER_LENGTH], h->BSIZE - Constants::TGZF_BLOCK_HEADER_LENGTH);
		if(!this->stream.good()){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Truncated file..." << std::endl;
			return false;
		}

		buffer.pointer = h->BSIZE;

		if(!this->gzip_controller.Inflate(buffer, output_buffer)){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Failed inflate!" << std::endl;
			return false;
		}

		if(this->output_buffer.size() == 0){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Empty data!" << std::endl;
			return false;
		}

		// Reset buffer
		this->buffer.reset();

		// Reset iterator position and size
		this->size = this->output_buffer.size() / sizeof(entry_type);

		// Validity check
		if(this->output_buffer.size() % sizeof(entry_type) != 0){
			std::cerr << Tomahawk::Helpers::timestamp("ERROR", "TWO") << "Data is corrupted!" << std::endl;
			return false;
		}
	}

	//std::cerr << "Loaded: " << this->output_buffer.size() << " pos " << this->stream.tellg() << std::endl;

	return true;
}

bool TomahawkOutputReader::nextVariant(const entry_type*& entry){
	if(this->position == this->size){
		if(!this->nextBlock())
			return false;
	}

	entry = (*this)[this->position];
	++this->position;

	return true;
}

bool TomahawkOutputReader::nextVariantLimited(const entry_type*& entry){
	if(this->position == this->size){
		return false;
	}

	entry = (*this)[this->position];
	++this->position;

	return true;
}


bool TomahawkOutputReader::summary(const std::string& input, const U32 bins){
	TWO::TomahawkOutputStatsContainer container(bins);

	// Natural output required parsing
	const entry_type* entry;
	while(this->nextVariant(entry))
		container += *entry;

	std::cerr << "R2\t" << container.R2.within.getTotal() << '\t' << container.R2.across.getTotal() << '\t' << container.R2.global.getTotal()  << std::endl;
	std::cerr << container.R2 << std::endl;
	std::cerr << "D\t" << container.D.within.getTotal() << '\t' << container.D.across.getTotal() << '\t' << container.D.global.getTotal()  << std::endl;
	std::cerr << container.D << std::endl;
	std::cerr << "Dprime\t" << container.Dprime.within.getTotal() << '\t' << container.Dprime.across.getTotal() << '\t' << container.Dprime.global.getTotal()  << std::endl;
	std::cerr << container.Dprime << std::endl;

	return true;
}

bool TomahawkOutputReader::index(const std::string& input){
	//if(!this->reader.setup(input))
	//	return false;

	const entry_type* entry;
	//if(!this->reader.nextEntry(entry))
	//	return false;

	//const entry_type* previous = entry;
	U32 currentAID = entry->AcontigID;
	U32 currentAPos = entry->Aposition;
	U32 currentBID = entry->BcontigID;
	U64 AIDSteps = 0;
	U64 APosSteps = 0;
	U64 BIDSteps = 0;

	double AposStepsR = 0;
	U64 outputEntries = 0;

	//while(this->reader.nextEntry(entry)){
	while(true){
		/*
		if(!(*previous < *entry)){
			std::cerr << "file is not sorted" << std::endl;
			std::cerr << previous->AcontigID << '\t' << previous->Aposition << '\t' << previous->BcontigID << '\t' << previous->Bposition << std::endl;
			std::cerr << entry->AcontigID << '\t' << entry->Aposition << '\t' << entry->BcontigID << '\t' << entry->Bposition << std::endl;
			return false;
		}
		std::swap(entry, previous);
		*/
		++AIDSteps;
		++APosSteps;
		++BIDSteps;
		AposStepsR += entry->R2;

		if(entry->BcontigID != currentBID || entry->Aposition != currentAPos || entry->AcontigID != currentAID){
			std::cerr << "switch: " << currentAID << ',' << currentAPos << ',' << currentBID << '\t' << entry->AcontigID << ',' << entry->Aposition << ',' << entry->BcontigID << '\t' << BIDSteps << std::endl;
			currentBID = entry->BcontigID;
			BIDSteps = 0;
			++outputEntries;
		}

		if(entry->Aposition != currentAPos || entry->AcontigID != currentAID){
			std::cout << currentAID << '\t' << currentAPos << '\t' << APosSteps << '\t' << AposStepsR/APosSteps << '\n';
			currentAPos = entry->Aposition;
			APosSteps = 0;
			AposStepsR = 0;
			++outputEntries;
		}

		if(entry->AcontigID != currentAID){
			std::cerr << "switch: " << currentAID << "->" << entry->AcontigID << '\t' << AIDSteps << std::endl;
			currentAID = entry->AcontigID;
			AIDSteps = 0;
			++outputEntries;
		}
	}
	std::cerr << "Index would have: " << outputEntries << " entries..." << std::endl;

	return true;
}

bool TomahawkOutputReader::setWriterType(const int type){
	if(type == 0)
		this->writer_output_type = WRITER_TYPE::binary;
	else if(type == 1)
		this->writer_output_type = WRITER_TYPE::natural;
	else {
		std::cerr << Tomahawk::Helpers::timestamp("ERROR","READER") << "Unknown writer type: " << type << std::endl;
		return false;
	}
	return true;
}


}
} /* namespace Tomahawk */
