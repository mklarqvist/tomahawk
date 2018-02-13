#include "index.h"

namespace Tomahawk{

Index::Index(){}
Index::~Index(){}

// Reading an index from a byte stream
Index::Index(const char* const data, const U32 l_data) :
	controller_(data),
	meta_container_(&data[sizeof(BYTE)], *reinterpret_cast<const size_type* const>(&data[sizeof(BYTE)])*TWK_INDEX_META_ENTRY_SIZE+sizeof(size_type)),
	container_(&data[sizeof(BYTE) + this->meta_container_.size() * TWK_INDEX_META_ENTRY_SIZE + sizeof(size_type)], l_data - (this->meta_container_.size() * TWK_INDEX_META_ENTRY_SIZE + sizeof(size_type) + sizeof(BYTE)))
{

}

bool Index::buildMetaIndex(const U32 n_contigs){
	if(this->getContainer().size() == 0)
		return false;

	if(this->isSorted() == false)
		return false;

	meta_entry_type reference_entry;
	reference_entry.index_begin = 0;
	reference_entry.index_end   = 1;
	reference_entry.min_position = this->getContainer()[0].min_position;
	reference_entry.max_position = this->getContainer()[0].max_position;
	reference_entry.n_variants = this->getContainer()[0].n_variants;
	reference_entry.uncompressed_size = this->getContainer()[0].uncompressed_size;
	U32 reference_contig = this->getContainer()[0].contigID;

	meta_entry_type* temp_entries = new meta_entry_type[n_contigs];

	for(U32 i = 1; i < this->getContainer().size(); ++i){
		if(this->getContainer()[i].contigID != reference_contig){
			if(this->getContainer()[i].contigID < reference_contig)
				continue;

			temp_entries[reference_contig] = reference_entry;
			reference_contig = this->getContainer()[i].contigID;
			reference_entry.index_begin = i;
			reference_entry.index_end = i + 1;
			reference_entry.min_position = this->getContainer()[i].min_position;
			reference_entry.max_position = this->getContainer()[i].max_position;
			reference_entry.n_variants = this->getContainer()[i].n_variants;
			reference_entry.uncompressed_size = this->getContainer()[i].uncompressed_size;

		} else {
			++reference_entry.index_end;
			reference_entry.max_position = this->getContainer()[i].max_position;
			reference_entry.n_variants += this->getContainer()[i].n_variants;
			reference_entry.uncompressed_size += this->getContainer()[i].uncompressed_size;
		}
	}
	temp_entries[reference_contig] = reference_entry;

	for(U32 i = 0; i < n_contigs; ++i){
		//std::cerr << temp_entries[i] << std::endl;
		this->meta_container_ += temp_entries[i];
	}
	delete [] temp_entries;

	return(true);
}

}
