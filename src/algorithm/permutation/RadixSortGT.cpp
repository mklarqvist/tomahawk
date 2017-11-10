#include "RadixSortGT.h"

namespace Tomahawk {
namespace Algorithm {

RadixSortGT::RadixSortGT() :
	n_samples(0),
	position(0),
	ppa(nullptr),
	GT_array(nullptr),
	bins(new U32*[9]),
	cumulative_AAC(0),
	cumulative_total(0),
	cost_ppa_conventional(0),
	cost_ppa_best(0),
	cost_ppa_byte(0),
	cost_ppa_u16(0),
	cost_ppa_u32(0),
	cost_ppa_u64(0),
	cost_rle_conventional(0),
	cost_rle_best(0),
	cost_rle_byte(0),
	cost_rle_u16(0),
	cost_rle_u32(0),
	cost_rle_u64(0)
{
	memset(&p_i, 0, sizeof(U32)*9);
}

RadixSortGT::RadixSortGT(const U64 n_samples) :
	n_samples(n_samples),
	position(0),
	ppa(new U32[this->n_samples]),
	GT_array(new BYTE[this->n_samples]),
	bins(new U32*[9]),
	cumulative_AAC(0),
	cumulative_total(0),
	cost_ppa_conventional(0),
	cost_ppa_best(0),
	cost_ppa_byte(0),
	cost_ppa_u16(0),
	cost_ppa_u32(0),
	cost_ppa_u64(0),
	cost_rle_conventional(0),
	cost_rle_best(0),
	cost_rle_byte(0),
	cost_rle_u16(0),
	cost_rle_u32(0),
	cost_rle_u64(0)
{
	memset(&p_i, 0, sizeof(U32)*9);
	for(U32 i = 0; i < 9; ++i){
		this->bins[i] = new U32[n_samples];
		memset(this->bins[i], 0, sizeof(U32)*n_samples);
	}

	memset(this->GT_array, 0, sizeof(BYTE)*n_samples);

	for(U32 i = 0; i < this->n_samples; ++i)
		this->ppa[i] = i;
}

RadixSortGT::~RadixSortGT(){
	delete [] this->ppa;
	delete [] this->GT_array;
	for(U32 i = 0; i < 9; ++i)
		delete [] this->bins[i];
}

void RadixSortGT::setSamples(const U64 n_samples){
	this->n_samples = n_samples;

	// Delete previous
	delete [] this->ppa;
	delete [] this->GT_array;

	// Set new
	this->ppa = new U32[this->n_samples];
	this->GT_array = new BYTE[this->n_samples];

	// Reset
	for(U32 i = 0; i < 9; ++i){
		this->bins[i] = new U32[n_samples];
		memset(this->bins[i], 0, sizeof(U32)*n_samples);
	}

	memset(this->GT_array, 0, sizeof(BYTE)*n_samples);

	for(U32 i = 0; i < this->n_samples; ++i)
		this->ppa[i] = i;
}

void RadixSortGT::reset(void){
	this->position = 0;
	memset(this->GT_array, 0, sizeof(BYTE)*n_samples);
	memset(&p_i, 0, sizeof(U32)*9);

	for(U32 i = 0; i < this->n_samples; ++i)
		this->ppa[i] = i;

	this->cumulative_AAC = 0;
	this->cumulative_total = 0;
}

bool RadixSortGT::build(const bcf_reader_type& reader){
	if(reader.size() == 0)
		return false;

	// Cycle over BCF entries
	for(U32 i = 0; i < reader.size(); ++i){
		// Has to be biallelic
		// otherwise skip
		if(!reader[i].isBiallelic())
			continue;

		//std::cerr << "update: " << i << std::endl;
		if(!this->update(reader[i]))
			continue;

	}

	// Return TRUE if the number of parsed
	// entries is > 0
	return(this->position > 0);
}

bool RadixSortGT::update(const bcf_entry_type& entry){
	// Check again because we might use it
	// iteratively at some point in time
	// i.e. not operating through the
	// build() function
	if(!entry.isBiallelic())
		return false;

	// Cycle over genotypes at this position
	// Ignore phasing at this stage
	//
	// Genotype encodings are thus:
	// 0/0 -> 0000b = 0 -> 0
	// 0/1 -> 0001b = 1 -> 3
	// 0/. -> 0010b = 2	-> 4
	// 1/0 -> 0100b = 4 -> 2
	// 1/1 -> 0101b = 5 -> 1
	// 1/. -> 0110b = 6 -> 5
	// ./0 -> 1000b = 8 -> 6
	// ./1 -> 1001b = 9 -> 7
	// ./. -> 1010b = 10 -> 8
	//
	// Update GT_array
	U32 alt = 0, ref = 0;
	U32 internal_pos = entry.p_genotypes;
	U32 k = 0;
	for(U32 i = 0; i < 2*this->n_samples; i += 2, ++k){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&entry.data[internal_pos++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&entry.data[internal_pos++]);
		const BYTE packed = (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) << 2) | BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1);
		this->GT_array[k] = packed;
		if(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) == 1) ++alt;
		else ++ref;
		if(BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1) == 1) ++alt;
		else ++ref;
	}
	if(alt <= ref) this->cumulative_AAC += alt;
	else this->cumulative_AAC += ref;
	this->cumulative_total += alt + ref;

	//if(alt < 50) return false;

	// Build PPA
	// 3^2 = 9 state radix sort over
	// states: alleles \in {00, 01, 11}
	// b entries in a TWK block B
	// This is equivalent to a radix sort
	// on the alphabet {0,1,...,8}
	U32 target_ID = 0;
	for(U32 j = 0; j < this->n_samples; ++j){
		// Determine correct bin
		switch(this->GT_array[this->ppa[j]]){
		case 0:  target_ID = 0; break;
		case 1:  target_ID = 3; break;
		case 2:  target_ID = 4; break;
		case 4:  target_ID = 2; break;
		case 5:  target_ID = 1; break;
		case 6:  target_ID = 5; break;
		case 8:  target_ID = 6; break;
		case 9:  target_ID = 7; break;
		case 10: target_ID = 8; break;
		default: std::cerr << "illegal" << std::endl; exit(1);
		}

		// Update bin i at position i with ppa[j]
		this->bins[target_ID][this->p_i[target_ID]] = this->ppa[j];
		++this->p_i[target_ID];
	} // end loop over individuals at position i

	// Update PPA data
	// Copy data in sorted order
	U32 cum_pos = 0;
	for(U32 i = 0; i < 9; ++i){
		// Copy data in bin i to current position
		memcpy(&this->ppa[cum_pos], this->bins[i], this->p_i[i]*sizeof(U32));

		// Update cumulative position and reset
		cum_pos += this->p_i[i];
		this->p_i[i] = 0;
	}
	// Make sure the cumulative position
	// equals the number of samples in the
	// dataset
	assert(cum_pos == this->n_samples);

	// Keep track of how many entries we've iterated over
	++this->position;

	return true;
}

void RadixSortGT::outputGT(const bcf_reader_type& reader){
	if(reader.size() == 0)
		return;

	for(U32 i = 0; i < reader.size(); ++i){
		if(!reader[i].isBiallelic())
			continue;

		U32 internal_pos = reader[i].p_genotypes;
		for(U32 k = 0; k < 2*this->n_samples; k += 2, ++k){
			const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&reader[i].data[internal_pos++]);
			const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&reader[i].data[internal_pos++]);
			const BYTE packed = (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) << 2) | BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1);
			this->GT_array[k] = packed;
		}

		for(U32 k = 0; k < this->n_samples; ++k){
			std::cerr << (int)this->GT_array[this->ppa[k]];
		}
		std::cerr << std::endl;
	}
}

bool RadixSortGT::assesRLECost(const bcf_reader_type& reader){
	if(reader.size() == 0)
		return false;

	for(U32 i = 0; i < reader.size(); ++i){
		if(!reader[i].isBiallelic())
			continue;

		// Load GT array
		U32 internal_pos = reader[i].p_genotypes;
		U32 j = 0;
		for(U32 k = 0; k < 2*this->n_samples; k += 2, ++j){
			const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&reader[i].data[internal_pos++]);
			const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&reader[i].data[internal_pos++]);
			const BYTE packed = (BCF::BCF_UNPACK_GENOTYPE(fmt_type_value2) << 2) | BCF::BCF_UNPACK_GENOTYPE(fmt_type_value1);
			this->GT_array[j] = packed;
		}

		this->assessRLEPPA(reader[i]);
		this->assessRLE(reader[i]);
	}

	this->cost_ppa_conventional += sizeof(U32)*this->n_samples;
	this->cost_ppa_best += sizeof(U32)*this->n_samples;
	this->cost_ppa_byte += sizeof(U32)*this->n_samples;
	this->cost_ppa_u16  += sizeof(U32)*this->n_samples;
	this->cost_ppa_u32  += sizeof(U32)*this->n_samples;
	this->cost_ppa_u64  += sizeof(U32)*this->n_samples;

	std::cerr << this->cost_ppa_conventional << '\t' << this->cost_ppa_best << '\t' << this->cost_ppa_byte << '\t' << this->cost_ppa_u16 << '\t' << this->cost_ppa_u32 << '\t' << this->cost_ppa_u64 << '\t' <<
				 this->cost_rle_conventional << '\t' << this->cost_rle_best << '\t' << this->cost_rle_byte << '\t' << this->cost_rle_u16 << '\t' << this->cost_ppa_u32 << '\t' << this->cost_rle_u64 << '\t' <<
				 double(this->cost_rle_best) / (this->cost_ppa_best) << '\t' << (signed long long)this->cost_rle_best - (this->cost_ppa_best) <<  '\t' <<
				 (double)this->cost_rle_conventional / this->cost_ppa_best << std::endl;


	return true;
}

U64 RadixSortGT::assessRLEPPA(const bcf_entry_type& entry){
	if(!entry.isBiallelic())
		return 0;

	U32 n_runsBYTE = 0;
	U32 n_runsU16 = 0;
	U32 n_runsU32 = 0;
	U32 n_runsU64 = 0;
	U32 n_runsNormal = 0;

	U32 current_BYTE = 1;
	U32 current_U16 = 1;
	U32 current_U32 = 1;
	U32 current_U64 = 1;
	BYTE ref = this->GT_array[this->ppa[0]];
	U32 longest_run = 1;

	for(U32 k = 1; k < this->n_samples; ++k){
		const BYTE internal_ref = this->GT_array[this->ppa[k]];
		if(ref != internal_ref){
			++n_runsBYTE;
			++n_runsU16;
			++n_runsU32;
			++n_runsU64;
			++n_runsNormal;
			if(current_BYTE > longest_run) longest_run = current_BYTE;
			if(current_U16 > longest_run) longest_run = current_U16;
			if(current_U32 > longest_run) longest_run = current_U32;
			if(current_U64 > longest_run) longest_run = current_U64;
			current_BYTE = 0;
			current_U16 = 0;
			current_U32 = 0;
			current_U64 = 0;
			ref = internal_ref;
		}
		if(current_BYTE == 63){
			current_BYTE = 0;
			++n_runsBYTE;
		}
		if(current_U16 == 16383){
			current_U16 = 0;
			++n_runsU16;
		}
		if(current_U32 == 1073741823){
			current_U32 = 0;
			++n_runsU32;
		}
		++current_BYTE;
		++current_U16;
		++current_U32;
		++current_U64;
	}
	++n_runsBYTE;
	++n_runsU16;
	++n_runsU32;
	++n_runsU64;
	++n_runsNormal;
	if(current_BYTE > longest_run) longest_run = current_BYTE;
	if(current_U16 > longest_run) longest_run = current_U16;
	if(current_U32 > longest_run) longest_run = current_U32;
	if(current_U64 > longest_run) longest_run = current_U64;

	BYTE conventional_cost = ceil((ceil(log2(double(longest_run))) + 2)/8);
	if(conventional_cost >= 3 && conventional_cost <= 4) conventional_cost = 4;
	else if(conventional_cost > 4) conventional_cost = 8;

	U64 smallest_cost = n_runsBYTE*sizeof(BYTE);
	if(n_runsU16*sizeof(U16) < smallest_cost) smallest_cost = n_runsU16*sizeof(U16);
	else if(n_runsU32*sizeof(U32) < smallest_cost) smallest_cost = n_runsU32*sizeof(U32);
	else if(n_runsU64*sizeof(U64) < smallest_cost) smallest_cost = n_runsU64*sizeof(U64);

	this->cost_ppa_conventional += conventional_cost * n_runsNormal;
	this->cost_ppa_best += smallest_cost;
	this->cost_ppa_byte += n_runsBYTE*sizeof(BYTE);
	this->cost_ppa_u16 += n_runsU16*sizeof(U16);
	this->cost_ppa_u32 += n_runsU32*sizeof(U32);
	this->cost_ppa_u64 += n_runsU64*sizeof(U64);

	return 0;
}

U64 RadixSortGT::assessRLE(const bcf_entry_type& entry){
	if(!entry.isBiallelic())
		return 0;

	U32 n_runsBYTE = 0;
	U32 n_runsU16 = 0;
	U32 n_runsU32 = 0;
	U32 n_runsU64 = 0;
	U32 n_runsNormal = 0;

	U32 current_BYTE = 1;
	U32 current_U16 = 1;
	U32 current_U32 = 1;
	U32 current_U64 = 1;
	BYTE ref = this->GT_array[0];
	U32 longest_run = 1;

	for(U32 k = 1; k < this->n_samples; ++k){
		const BYTE internal_ref = this->GT_array[k];
		if(internal_ref != ref){
			++n_runsBYTE;
			++n_runsU16;
			++n_runsU32;
			++n_runsU64;
			++n_runsNormal;
			if(current_BYTE > longest_run) longest_run = current_BYTE;
			if(current_U16 > longest_run) longest_run = current_U16;
			if(current_U32 > longest_run) longest_run = current_U32;
			if(current_U64 > longest_run) longest_run = current_U64;
			current_BYTE = 0;
			current_U16 = 0;
			current_U32 = 0;
			current_U64 = 0;
			ref = internal_ref;
		}
		if(current_BYTE == 63){
			current_BYTE = 0;
			++n_runsBYTE;
		}
		if(current_U16 == 16383){
			current_U16 = 0;
			++n_runsU16;
		}
		if(current_U32 == 1073741823){
			current_U32 = 0;
			++n_runsU32;
		}
		++current_BYTE;
		++current_U16;
		++current_U32;
		++current_U64;
	}
	++n_runsBYTE;
	++n_runsU16;
	++n_runsU32;
	++n_runsU64;
	++n_runsNormal;
	if(current_BYTE > longest_run) longest_run = current_BYTE;
	if(current_U16 > longest_run) longest_run = current_U16;
	if(current_U32 > longest_run) longest_run = current_U32;
	if(current_U64 > longest_run) longest_run = current_U64;

	BYTE conventional_cost = ceil((ceil(log2(double(longest_run))) + 2)/8);
	if(conventional_cost >= 3 && conventional_cost <= 4) conventional_cost = 4;
	else if(conventional_cost > 4) conventional_cost = 8;

	U64 smallest_cost = n_runsBYTE*sizeof(BYTE);
	if(n_runsU16*sizeof(U16) < smallest_cost) smallest_cost = n_runsU16*sizeof(U16);
	else if(n_runsU32*sizeof(U32) < smallest_cost) smallest_cost = n_runsU32*sizeof(U32);
	else if(n_runsU64*sizeof(U64) < smallest_cost) smallest_cost = n_runsU64*sizeof(U64);

	this->cost_rle_conventional += conventional_cost * n_runsNormal;
	this->cost_rle_best += smallest_cost;
	this->cost_rle_byte += n_runsBYTE*sizeof(BYTE);
	this->cost_rle_u16 += n_runsU16*sizeof(U16);
	this->cost_rle_u32 += n_runsU32*sizeof(U32);
	this->cost_rle_u64 += n_runsU64*sizeof(U64);

	return 0;
}


} /* namespace IO */
} /* namespace Tomahawk */
