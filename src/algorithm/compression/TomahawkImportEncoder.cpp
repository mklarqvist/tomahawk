#include "TomahawkImportEncoder.h"

namespace Tomahawk{
namespace Algorithm{

bool TomahawkImportEncoder::Encode(const bcf_type& line, meta_base_type& meta_base, buffer_type& runs, buffer_type& simple, U64& n_runs){
	if(line.body->n_allele + 1 >= 32768){
		std::cerr << Helpers::timestamp("ERROR", "ENCODER") <<
					 "Illegal number of alleles (" << line.body->n_allele + 1 << "). "
					 "Format is limited to 32768..." << std::endl;
		return false;
	}

	// Update basic values for a meta entry
	meta_base.position = (U32)line.body->POS + 1; // Base-1
	meta_base.ref_alt = line.ref_alt;
	meta_base.controller.biallelic = true;
	meta_base.controller.simple = line.isSimple();

	// Assess cost and encode
	rle_helper_type cost;
	if(line.body->n_allele == 2){
		cost = this->assessRLEBiallelic(line);
		meta_base.virtual_offset_gt = runs.pointer; // absolute position at start of stream
		meta_base.controller.rle = true;
		meta_base.controller.mixed_phasing = cost.mixedPhasing;
		meta_base.controller.anyMissing = cost.hasMissing;

		//std::cerr << line.body->POS+1 << ":" << (int)cost.word_width << ',' << cost.n_runs << ',' << (int)cost.hasMissing << ',' << (int)cost.mixedPhasing << std::endl;

		switch(cost.word_width){
		case 1:
			this->EncodeRLE<BYTE>(line, runs, n_runs, cost.hasMissing, cost.mixedPhasing);
			meta_base.controller.rle_type = 0;
			break;
		case 2:
			this->EncodeRLE<U16>(line, runs, n_runs, cost.hasMissing, cost.mixedPhasing);
			meta_base.controller.rle_type = 1;
			break;
		case 4:
			this->EncodeRLE<U32>(line, runs, n_runs, cost.hasMissing, cost.mixedPhasing);
			meta_base.controller.rle_type = 2;
			break;
		case 8:
			this->EncodeRLE<U64>(line, runs, n_runs, cost.hasMissing, cost.mixedPhasing);
			meta_base.controller.rle_type = 3;
			break;
		default:
			std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
			return false;
		}

		meta_base.AF = float(this->helper.countsAlleles[1]) / (this->helper.countsAlleles[0] + this->helper.countsAlleles[1]);
		// See meta for more information
		meta_base.MGF = this->helper.MGF;
		meta_base.HWE_P = this->helper.HWE_P;

		// Reset and recycle helper
		this->helper.reset();

		return true;
	}
	else {
		cost = this->assessRLEnAllelic(line);
		U32 costBCFStyle = this->n_samples;
		if(line.body->n_allele + 1 < 8)          costBCFStyle *= 1;
		else if(line.body->n_allele + 1 < 128)   costBCFStyle *= 2;
		else if(line.body->n_allele + 1 < 32768) costBCFStyle *= 4;

		meta_base.virtual_offset_gt = simple.pointer; // absolute position at start of stream
		meta_base.controller.biallelic = false;
		meta_base.AF = 0; // AF needs to be looked up in cold store
		meta_base.controller.mixed_phasing = cost.mixedPhasing;
		meta_base.controller.anyMissing = cost.hasMissing;

		// See meta for more information
		meta_base.MGF = 0;
		meta_base.HWE_P = 1;

		//std::cerr << (int)cost.word_width << '\t' << cost.n_runs << '\t' << cost.word_width*cost.n_runs << '\t' << costBCFStyle << '\t' << costBCFStyle/double(cost.word_width*cost.n_runs) << std::endl;

		// RLE is cheaper
		if(cost.word_width*cost.n_runs < costBCFStyle){
			meta_base.controller.rle = true;
			meta_base.controller.rle_type = 1; // Todo: fix

			switch(cost.word_width){
			case 1:
				this->EncodeRLESimple<BYTE>(line, runs, n_runs);
				meta_base.controller.rle_type = 0;
				break;
			case 2:
				this->EncodeRLESimple<U16>(line, runs, n_runs);
				meta_base.controller.rle_type = 1;
				break;
			case 4:
				this->EncodeRLESimple<U32>(line, runs, n_runs);
				meta_base.controller.rle_type = 2;
				break;
			case 8:
				this->EncodeRLESimple<U64>(line, runs, n_runs);
				meta_base.controller.rle_type = 3;
				break;
			default:
				std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
				return false;
			}

			// Reset and recycle helper
			this->helper.reset();
			return true;
		}
		// BCF style is cheaper
		else {
			meta_base.controller.rle = false;
			meta_base.controller.rle_type = 0;

			if(line.body->n_allele + 1 < 8)          this->EncodeSingle<BYTE>(line, simple, n_runs);
			else if(line.body->n_allele + 1 < 128)   this->EncodeSingle<U16> (line, simple, n_runs);
			else if(line.body->n_allele + 1 < 32768) this->EncodeSingle<U32> (line, simple, n_runs);
			else {
				std::cerr << Helpers::timestamp("ERROR", "ENCODER") <<
							 "Illegal number of alleles (" << line.body->n_allele + 1 << "). "
							 "Format is limited to 32768..." << std::endl;
				return false;
			}

			// Reset and recycle helper
			this->helper.reset();
			return true;
		}
	}
	return false;
}

const TomahawkImportEncoder::rle_helper_type TomahawkImportEncoder::assessRLEBiallelic(const bcf_type& line){
	// Assess RLE cost
	U32 internal_pos_rle = line.p_genotypes;
	const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
	const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
	U32 n_runs = 1;
	U32 largest = 0;
	U32 run_length = 1;
	bool mixedPhase = false;
	bool anyMissing = false;
	const BYTE firstPhase = (fmt_type_value2 & 1);

	U32 ref = PACK_RLE_BIALLELIC(fmt_type_value2, fmt_type_value1, 2, 1);

	for(U32 i = 2; i < this->n_samples * 2; i += 2){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
		U32 internal = PACK_RLE_BIALLELIC(fmt_type_value2, fmt_type_value1, 2, 1);

		// Any data is missing?
		if((fmt_type_value1 >> 1) == 0 || (fmt_type_value2 >> 1) == 0) anyMissing = true;
		// Is there mixed phasing?
		if((fmt_type_value2 & 1) != firstPhase) mixedPhase = true;

		// Extend or break run
		if(ref != internal){
			ref = internal;
			if(run_length > largest) largest = run_length;
			++n_runs;
			run_length = 0;
		}
		++run_length;
	}
	++n_runs;
	if(run_length > largest) largest = run_length;

	// Calculate word size
	// Allele width can accidentally be 0 if log2(1)
	BYTE allele_width = 2*ceil(log2(1 + anyMissing));
	if(allele_width == 0) allele_width = 1 + 1;
	U32 word_width = ceil((ceil(log2(largest + 1)) + allele_width + mixedPhase) / 8);

	if(word_width == 3) word_width = 4;
	if(word_width > 4)  word_width = 8;
	return(rle_helper_type(word_width, n_runs, mixedPhase, anyMissing));
}

const TomahawkImportEncoder::rle_helper_type TomahawkImportEncoder::assessRLEnAllelic(const bcf_type& line){
	// Assess RLE cost
	const BYTE shift_size = ceil(log2(line.body->n_allele + 1));
	U32 internal_pos_rle = line.p_genotypes;
	const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
	const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
	U32 n_runs = 1;
	U32 largest = 0;
	U32 run_length = 1;
	U32 ref = ((fmt_type_value2 >> 1) << (shift_size + 1)) |
			  ((fmt_type_value1 >> 1) << 1) |
			  (fmt_type_value2 & 1);

	for(U32 i = 2; i < this->n_samples * 2; i += 2){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
		U32 internal = ((fmt_type_value2 >> 1) << (shift_size + 1)) |
					   ((fmt_type_value1 >> 1) << 1) |
					   (fmt_type_value2 & 1);

		if(ref != internal){
			ref = internal;
			if(run_length > largest) largest = run_length;
			++n_runs;
			run_length = 0;
		}
		++run_length;
	}
	++n_runs;
	if(run_length > largest) largest = run_length;

	// Calculate word size
	// Allele width can accidentally be 0 if log2(1)
	U32 word_width = ceil(
			         (ceil(log2(largest + 1)) +
			         2*ceil(log2(line.body->n_allele + 1)) + 1)
			         / 8);
	if(word_width == 3) word_width = 4;
	if(word_width > 4)  word_width = 8;
	return(rle_helper_type(word_width, n_runs, true, true));
}

}
}
