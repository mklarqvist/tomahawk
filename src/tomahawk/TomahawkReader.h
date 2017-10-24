#ifndef TOMAHAWK_TOMAHAWKREADER_H_
#define TOMAHAWK_TOMAHAWKREADER_H_

#include <string>
#include <fstream>
#include <algorithm>
#include <bitset>

#include "../support/MagicConstants.h"
#include "../io/compression/TGZFController.h"
#include "../io/compression/GZFConstants.h"
#include "../interface/Timer.h"
#include "../interface/ProgressBar.h"
#include "../algorithm/Balancer.h"
#include "TomahawkCalcParameters.h"
#include "base/TomahawkEntryMeta.h"
#include "TomahawkCalculateSlave.h"

namespace Tomahawk {
namespace Support{

// Data structure used together with the Occ
// function. Keeps diploid genotype lookups in
// a simple data structure
struct GroupGenotypes{
	GroupGenotypes(void) : count(0){
		memset(&genotypes[0], 0, sizeof(U64)*16);
	}

	~GroupGenotypes(){}

	void reset(void){
		if(count == 0)
			return;

		memset(&genotypes[0], 0, sizeof(U64)*16);
		count = 0;
	}

	void add(const BYTE& genotype, const U64& length){
		this->genotypes[genotype] += length;
		count += length;
	}

	inline U64& operator[](const U32& p){ return(this->genotypes[p]); }
	inline const U64& operator[](const U32& p) const{ return(this->genotypes[p]); }

	U64 count;
	U64 genotypes[16];
};

}

// TomahawkReader class simply reads compressed data from disk
class TomahawkReader {
	typedef TomahawkCalcParameters parameter_type;
	typedef Totempole::TotempoleEntry totempole_entry;
	typedef IO::TGZFController tgzf_controller_type;
	typedef Totempole::TotempoleReader totempole_type;
	typedef IO::BasicBuffer buffer_type;
	typedef IO::GenericWriterInterace writer_interface;
	typedef std::vector<U64> occ_vector;
	typedef std::vector<occ_vector> occ_matrix;
	typedef Tomahawk::Hash::HashTable<std::string, S32> hash_table;

public:
	// Used to keep track of char pointer offsets in buffer
	// and what Totempole entry is associated with that position
	// given the data that is loaded
	// This is equivalent to Totempole virtual offsets given
	// the data loaded in memory
	struct DataOffsetPair{
		DataOffsetPair(const char* data, const totempole_entry& entry) : entry(entry), data(data){}
		~DataOffsetPair(){}

		const totempole_entry& entry;
		const char* data;
	};

	struct GroupPair{
		GroupPair(const std::string& name) : n_entries(1), name(name){}
		~GroupPair(){}

		void operator++(void){ ++this->n_entries; }
		void operator--(void){ --this->n_entries; }
		void operator+=(const U32 p){ this->n_entries += p; }
		void operator-=(const U32 p){ this->n_entries -= p; }

		U32 n_entries;
		std::string name;
	};

private:
	typedef std::vector<DataOffsetPair> offset_vector;

public:
	TomahawkReader();
	virtual ~TomahawkReader();

	bool Open(const std::string input);
	bool loadGroups(const std::string& file);

	// Reader functions
	bool nextBlock(const bool clear = true);
	bool getBlocks(void);
	bool getBlocks(std::vector<U32>& blocks);
	bool getBlocks(std::vector< std::pair<U32, U32> >& blocks);
	bool getBlock(const U32 blockID);

	// Output functions
	bool outputBlocks(std::vector<U32>& blocks);
	bool outputBlocks();

	inline const BYTE& getBitWidth(void) const{ return(this->bit_width); }
	inline Totempole::TotempoleReader& getTotempole(void){ return(this->totempole); }
	inline const DataOffsetPair& getOffsetPair(const U32 p) const{ return(this->blockDataOffsets[p]); }
	inline const size_t DataOffsetSize(void) const{ return(this->blockDataOffsets.size()); }
	inline void setDropGenotypes(const bool yes){ this->dropGenotypes = yes; }
	inline void setShowHeader(const bool yes){ this->showHeader = yes; }

private:
	void DetermineBitWidth(void);
	bool Validate(void);
	bool ValidateHeader(std::ifstream& in) const;
	template <class T> bool __calculateTajimaD(const U32 bin_size);
	template <class T> bool __calculateFST(void);
	template <class T> bool __calculateSFS(void);
	template <class T> bool __outputBlock(const U32& id);
	template <class T> bool __outputBlockGrouped(const U32& id);

protected:
	U64 samples;    // has to match header
	float version;  // has to match header
	U64 filesize;   // filesize
	BYTE bit_width; // bit width
	bool dropGenotypes; // drop genotypes in view mode
	bool showHeader;// flag to output header or not

	U32 currentBlockID; // for iterator

	std::ifstream stream; // reader stream

	totempole_type totempole; // totempole reader
	buffer_type buffer; // input buffer
	buffer_type data; // inflate buffer
	buffer_type outputBuffer; // output buffer
	tgzf_controller_type tgzf_controller; // tgzf controller
	offset_vector blockDataOffsets; // internal virtual offsets into buffer
	writer_interface* writer; // writer interface

	// groups
	occ_matrix Occ;
	std::vector<GroupPair> groups;
	hash_table* group_htable;
};

template <class T>
bool TomahawkReader::__outputBlock(const U32& id){
	TomahawkIterator<T> tomahawk_controller(this->data.data, this->totempole[id]);

	// For each variant in Tomahawk block
	for(U32 j = 0; j < tomahawk_controller.support->n_variantsRLE; ++j){
		tomahawk_controller.WriteVariant(this->totempole, this->outputBuffer, this->dropGenotypes);

		// Next variant
		++tomahawk_controller;

		// Keep flushing regularly
		// arbitrary threshold at 65536 bytes
		if(this->outputBuffer.size() > 65536){
			//this->writer->write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
			std::cout.write(&this->outputBuffer.data[0], this->outputBuffer.pointer);
			this->outputBuffer.reset();
		}
	}

	// Flush last
	//this->writer->write(&this->outputBuffer_.data[0], this->outputBuffer_.pointer);
	std::cout.write(&this->outputBuffer.data[0], this->outputBuffer.pointer);

	// Reset buffers
	this->outputBuffer.reset(); // reset
	this->data.reset(); // reset

	return true;
}

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKREADER_H_ */
