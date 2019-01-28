#include "ld_engine.h"
#include "fisher_math.h"

namespace tomahawk {

twk_ld_simd::twk_ld_simd(void) :
	counters((uint64_t*)_mm_malloc(sizeof(uint64_t)*16, 16))
{
	memset(this->counters, 0, sizeof(uint64_t)*16);
}

twk_ld_simd::~twk_ld_simd(){
	_mm_free(this->counters);
}

twk_ld_count::twk_ld_count(): totalHaplotypeCounts(0)
{
	// Initialize counters to 0. This is generally not necessary as each
	// function computing LD clears these. However, it is good practice.
	memset(alleleCounts,  171, sizeof(uint64_t)*171);
	memset(haplotypeCounts, 4, sizeof(uint64_t)*4);
}
twk_ld_count::~twk_ld_count(){}

void twk_ld_count::ResetPhased(void){
	memset(alleleCounts,    0, sizeof(uint64_t)*16);
	memset(haplotypeCounts, 0, sizeof(uint64_t)*4);
}

void twk_ld_count::ResetUnphased(void){
	memset(alleleCounts,    0, sizeof(uint64_t)*171);
	memset(haplotypeCounts, 0, sizeof(uint64_t)*4);
}

//
twk_ld_engine::twk_ld_engine() :
	n_samples(0), n_out(0), n_lim(10000), n_out_tick(250),
	byte_width(0), byte_aligned_end(0), vector_cycles(0),
	phased_unbalanced_adjustment(0), unphased_unbalanced_adjustment(0), t_out(0),
	mask_placeholder(nullptr),
	index(nullptr), writer(nullptr), progress(nullptr), list_out(nullptr)
{
	memset(n_method, 0, sizeof(uint64_t)*10);
}

twk_ld_engine::~twk_ld_engine(){ delete[] list_out; aligned_free(mask_placeholder); }

/**<
 * Set the number of samples in the target file. This function is mandatory
 * to call prior to computing as it sets essential paramters such as the number
 * of samples, the byte alignment of vectorized bit-vectors and the adjustment
 * given the unused overhangs.
 * @param samples Number of samples in the target file.
 */
void twk_ld_engine::SetSamples(const uint32_t samples){
	n_samples  = samples;

	byte_width       = std::ceil(2.0f*samples/64);
	vector_cycles    = 2*samples/SIMD_WIDTH; // integer division
	byte_aligned_end = vector_cycles * (SIMD_WIDTH / 64);
	phased_unbalanced_adjustment   = (byte_width*64 - 2*samples) / 2;
	unphased_unbalanced_adjustment = (byte_width*64 - 2*samples) / 2;

	//std::cerr << byte_width << "," << vector_cycles << "," << byte_aligned_end << "," << phased_unbalanced_adjustment << "," << unphased_unbalanced_adjustment << std::endl;
	//exit(1);

	uint32_t n = ceil((double)(n_samples*2)/64);
	n += (n*64) % SIMD_WIDTH; // must be divisible by 128-bit register
	mask_placeholder = reinterpret_cast<uint64_t*>(aligned_malloc(n*sizeof(uint64_t), SIMD_ALIGNMENT));
	memset(mask_placeholder, 0, n*sizeof(uint64_t));
}

void twk_ld_engine::SetBlocksize(const uint32_t s){
	assert(s % 2 == 0);
	blk_f.clear(); blk_r.clear();
	blk_f.resize(s + 100);
	blk_r.resize(s + 100);
	obuf.resize(s * sizeof(twk1_two_t));
	ibuf.resize(s * sizeof(twk1_two_t));

	n_out = 0; n_lim = s;
}

bool twk_ld_engine::PhasedList(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
	helper.ResetPhased();
	const twk_igt_list& ref = b1.list[p1];
	const twk_igt_list& tgt = b2.list[p2];
	const uint32_t n_cycles = ref.l_list < tgt.l_list ? ref.l_list : tgt.l_list;

	// Debug timings
	#if TWK_SLAVE_DEBUG_MODE == 1
		typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
		auto t0 = std::chrono::high_resolution_clock::now();
	#endif

	if(n_cycles == 0) helper.alleleCounts[TWK_LD_ALTALT] = 0;
	else {
		if(ref.l_list >= tgt.l_list){
			if(tgt.list[tgt.l_list - 1] < ref.list[0]) {
				helper.alleleCounts[TWK_LD_ALTALT] = 0;
				//std::cerr << "no overlap=" << tgt.list[tgt.l_list - 1] << "<" << ref.list[0] << " len=" << tgt.l_list << "<" << ref.l_list << std::endl;
			}
			else
				for(uint32_t i = 0; i < n_cycles; ++i) helper.alleleCounts[TWK_LD_ALTALT] += ref.get(tgt.list[i]);
		} else {
			if(ref.list[ref.l_list - 1] < tgt.list[0]) {
				helper.alleleCounts[TWK_LD_ALTALT] = 0;
				//std::cerr << "no overlap=" << ref.list[ref.l_list - 1] << "<" << tgt.list[0] << " len=" << ref.l_list << "<" << tgt.l_list << std::endl;
			} else
				for(uint32_t i = 0; i < n_cycles; ++i) helper.alleleCounts[TWK_LD_ALTALT] += tgt.get(ref.list[i]);
		}
	}

	helper.alleleCounts[TWK_LD_ALTREF] = ref.l_list - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFALT] = tgt.l_list - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFREF] = 2*n_samples - ((ref.l_list + tgt.l_list) - helper.alleleCounts[TWK_LD_ALTALT]);
	++n_method[0];

#if TWK_SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[n_cycles] += ticks_per_iter.count();
	++perf->freq[n_cycles];
#endif

#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "list=" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
#endif

#if TWK_SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::UnphasedList(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
	helper.ResetUnphased();

	const twk_igt_list& refBV = b1.list[p1];
	const twk_igt_list& tgtBV = b2.list[p2];

	// Step 1: count (alt,alt),(alt,alt) tuples
	// Step 2: count (het,het) tuples
	// Step 3: infer remainder

	// g(A) = popcount(((A >> 1) & A) & low_mask)
	// g(A) = popcnt128(helper.alleleCounts[TWK_LD_ALTALT], _mm_and_si128(_mm_and_si128(_mm_slli_epi64(PHASED_ALTALT(vectorA,vectorB), 1), PHASED_ALTALT(vectorA,vectorB)), maskUnphasedLow))
	// x = PHASED_ALTALT(vectorA,vectorB)

	// Compute joint alts and alt refs
	if(refBV.r_aa.size() < tgtBV.r_aa.size()){
		for(uint32_t i = 0; i < refBV.r_pos.size(); ++i){
			const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)&refBV.bv[refBV.r_aa[i]]);
			const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)&tgtBV.bv[refBV.r_aa[i]]);
			popcnt128(helper.alleleCounts[85], _mm_and_si128(_mm_and_si128(_mm_slli_epi64(PHASED_ALTALT(vectorA,vectorB), 1), PHASED_ALTALT(vectorA,vectorB)), maskUnphasedLow));
		}
	} else {
		for(uint32_t i = 0; i < tgtBV.r_aa.size(); ++i){
			const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)&refBV.bv[tgtBV.r_aa[i]]);
			const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)&tgtBV.bv[tgtBV.r_aa[i]]);
			popcnt128(helper.alleleCounts[85], _mm_and_si128(_mm_and_si128(_mm_slli_epi64(PHASED_ALTALT(vectorA,vectorB), 1), PHASED_ALTALT(vectorA,vectorB)), maskUnphasedLow));
		}
	}

	// Compute hets any and joint hets
	if(refBV.r_het.size() < tgtBV.r_het.size()){
		for(uint32_t i = 0; i < refBV.r_het.size(); ++i){
			const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)&refBV.bv[refBV.r_het[i]]);
			const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)&tgtBV.bv[refBV.r_het[i]]);
			popcnt128(helper.alleleCounts[TWK_LD_ALTALT], _mm_and_si128(vectorA,vectorB));
		}
	} else {
		for(uint32_t i = 0; i < tgtBV.r_het.size(); ++i){
			const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)&refBV.bv[tgtBV.r_het[i]]);
			const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)&tgtBV.bv[tgtBV.r_het[i]]);
			popcnt128(helper.alleleCounts[TWK_LD_ALTALT], _mm_and_si128(vectorA,vectorB));
		}
	}


	return false;
}

bool twk_ld_engine::PhasedListVector(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
	helper.ResetPhased();

	const twk_igt_list& refBV = b1.list[p1];
	const twk_igt_list& tgtBV = b2.list[p2];

	// Debug timings
#if TWK_SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	const VECTOR_TYPE binmapA = _mm_loadu_si128((const VECTOR_TYPE*)refBV.bin_bitmap);
	const VECTOR_TYPE binmapB = _mm_loadu_si128((const VECTOR_TYPE*)tgtBV.bin_bitmap);
	const VECTOR_TYPE binmap  = _mm_and_si128(binmapA,binmapB);
	//const uint64_t b[2] = {_mm_extract_epi64(binmap, 0), _mm_extract_epi64(binmap, 1)};

	// Check bins of bins if anything overlaps
	uint64_t a = 0;
	popcnt128(a, binmap);

	if(a == 0){
#if TWK_SLAVE_DEBUG_MODE == 1
		//const uint32_t n_list = refBV.l_list < tgtBV.l_list ? refBV.l_list : tgtBV.l_list;
		const uint32_t n_list = b1.list[p1].d.n < b2.list[p2].d.n ? b1.list[p1].d.n : b2.list[p2].d.n;
		auto t1 = std::chrono::high_resolution_clock::now();
		auto ticks_per_iter = Cycle(t1-t0);
		perf->cycles[n_list] += ticks_per_iter.count();
		++perf->freq[n_list];
#endif
		helper.alleleCounts[TWK_LD_ALTALT] = 0;
		helper.alleleCounts[TWK_LD_ALTREF] = b1.list[p1].l_list - helper.alleleCounts[TWK_LD_ALTALT];
		helper.alleleCounts[TWK_LD_REFALT] = b2.list[p2].l_list - helper.alleleCounts[TWK_LD_ALTALT];
		helper.alleleCounts[TWK_LD_REFREF] = 2*n_samples - ((b1.list[p1].l_list + b2.list[p2].l_list) - helper.alleleCounts[TWK_LD_ALTALT]);
#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "listBV-exit=" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
#endif
		++n_method[0];
#if TWK_SLAVE_DEBUG_MODE != 1
		return(PhasedMath(b1,p1,b2,p2));
#else
		return(true);
#endif
	}

	if(refBV.r_pos.size() < tgtBV.r_pos.size()){
		for(uint32_t i = 0; i < refBV.r_pos.size(); ++i){
			const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)&refBV.bv[refBV.r_pos[i]]);
			const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)&tgtBV.bv[refBV.r_pos[i]]);
			popcnt128(helper.alleleCounts[TWK_LD_ALTALT], _mm_and_si128(vectorA,vectorB));
		}
	} else {
		for(uint32_t i = 0; i < tgtBV.r_pos.size(); ++i){
			const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)&refBV.bv[tgtBV.r_pos[i]]);
			const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)&tgtBV.bv[tgtBV.r_pos[i]]);
			popcnt128(helper.alleleCounts[TWK_LD_ALTALT], _mm_and_si128(vectorA,vectorB));
		}
	}

	helper.alleleCounts[TWK_LD_ALTREF] = refBV.l_list - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFALT] = tgtBV.l_list - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFREF] = 2*n_samples - ((refBV.l_list + tgtBV.l_list) - helper.alleleCounts[TWK_LD_ALTALT]);
	++n_method[0];

#if TWK_SLAVE_DEBUG_MODE == 1
	//const uint32_t n_list = refBV.l_list < tgtBV.l_list ? refBV.l_list : tgtBV.l_list;
	const uint32_t n_list = b1.list[p1].d.n < b2.list[p2].d.n ? b1.list[p1].d.n : b2.list[p2].d.n;
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[n_list] += ticks_per_iter.count();
	++perf->freq[n_list];
#endif

#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "listBV=" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
#endif

#if TWK_SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

/*
bool twk_ld_engine::PhasedListSpecial(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
	helper.ResetPhased();
	//const twk_igt_list::ilist_cont& ref = b1.list[p1].d;
	//const twk_igt_list::ilist_cont& tgt = b2.list[p2].d;
	const twk_igt_list& refBV = b1.list[p1];
	const twk_igt_list& tgtBV = b2.list[p2];

	// Debug timings
#if TWK_SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	const VECTOR_TYPE binmapA = _mm_load_si128((const VECTOR_TYPE*)refBV.bin_bitmap);
	const VECTOR_TYPE binmapB = _mm_load_si128((const VECTOR_TYPE*)tgtBV.bin_bitmap);

	// Check bins of bins if anything overlaps
	uint64_t a = 0;
	popcnt128(a, _mm_and_si128(binmapA,binmapB));
	//std::cerr << "and=" << a << " " << std::bitset<64>(refBV.bin_bitmap[0] & tgtBV.bin_bitmap[0])  << " " << std::bitset<64>(refBV.bin_bitmap[1] & tgtBV.bin_bitmap[1])<< std::endl;

	if(a == 0){
#if TWK_SLAVE_DEBUG_MODE == 1
		auto t1 = std::chrono::high_resolution_clock::now();
		auto ticks_per_iter = Cycle(t1-t0);
		//const uint32_t n_list = refBV.l_list < tgtBV.l_list ? refBV.l_list : tgtBV.l_list;
		const uint32_t n_list = refBV.x_bins < tgtBV.x_bins ? refBV.x_bins : tgtBV.x_bins;
		perf->cycles[n_list] += ticks_per_iter.count();
		++perf->freq[n_list];
#endif
		helper.alleleCounts[TWK_LD_ALTALT] = 0;
		helper.alleleCounts[TWK_LD_ALTREF] = b1.list[p1].l_list - helper.alleleCounts[TWK_LD_ALTALT];
		helper.alleleCounts[TWK_LD_REFALT] = b2.list[p2].l_list - helper.alleleCounts[TWK_LD_ALTALT];
		helper.alleleCounts[TWK_LD_REFREF] = 2*n_samples - ((b1.list[p1].l_list + b2.list[p2].l_list) - helper.alleleCounts[TWK_LD_ALTALT]);
#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "listS-exit=" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
#endif
		++n_method[0];
#if TWK_SLAVE_DEBUG_MODE != 1
		return(PhasedMath(b1,p1,b2,p2));
#else
		return(true);
#endif
	}


	// Cycle over buckets of registers
	if(refBV.x.pos.size() < tgtBV.x.pos.size()){
		//std::cerr << "pos=" << refBV.x.pos.size() << std::endl;
		for(int k = 0; k < refBV.x.pos.size(); ++k){
			const uint32_t reg = refBV.x.pos[k];
			//std::cerr << "curpos=" << reg << std::endl;
			if(!(refBV.x.data[reg].s && tgtBV.x.data[reg].s)) continue;
				//uint32_t i = 0, j = 0;
				//uint32_t ii = 0, jj = 0;

				const twk_igt_list::ilist_cont& ref = refBV.x.data[reg];
				const twk_igt_list::ilist_cont& tgt = tgtBV.x.data[reg];
				const std::vector<uint32_t>& posA = ref.pos;
				const std::vector<uint32_t>& posB = tgt.pos;

				// Lookup tgt register in ref register
				if(posA.size() < posB.size()){
					for(int i = 0; i < posA.size(); ++i){
						const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)ref.data[posA[i]].vals);
						const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)tgt.data[posA[i]].vals);
						popcnt128(helper.alleleCounts[TWK_LD_ALTALT], _mm_and_si128(vectorA,vectorB));
					}
				} else {
					for(int i = 0; i < posB.size(); ++i){
						const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)ref.data[posB[i]].vals);
						const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)tgt.data[posB[i]].vals);
						popcnt128(helper.alleleCounts[TWK_LD_ALTALT], _mm_and_si128(vectorA,vectorB));
					}
				}

				/*while(true){
					// select
					const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)refBV.x.data[reg].data[i].vals);
					const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)tgtBV.x.data[reg].data[j].vals);

					popcnt128(helper.alleleCounts[TWK_LD_ALTALT],
							_mm_and_si128(
							  _mm_and_si128(vectorA,vectorB),
											_mm_set1_epi32(mask_lookup[refBV.x.data[reg].data[i].bin == tgtBV.x.data[reg].data[j].bin])));

					//std::cerr << "Diff@" << k << " -> " << refBV.x.data[k].data[i].bin << " and " << tgtBV.x.data[k].data[j].bin << std::endl;
					ii = i; jj = j;
					i += (refBV.x.data[reg].data[ii].bin == tgtBV.x.data[reg].data[jj].bin); // A == B
					j += (refBV.x.data[reg].data[ii].bin == tgtBV.x.data[reg].data[jj].bin); // A == B
					j += refBV.x.data[reg].data[ii].bin > tgtBV.x.data[reg].data[jj].bin; // A > B
					i += refBV.x.data[reg].data[ii].bin < tgtBV.x.data[reg].data[jj].bin; // A < B

					//std::cerr << i << "/" << ref.n << " and " << j << "/" << tgt.n << std::endl;
					if(i == refBV.x.data[reg].n || j == tgtBV.x.data[reg].n) break;
				}*\/
				//std::cerr << "overlap=" << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
			//}
			//else {
				//std::cerr << "no overlap = " << refBV.x.data[k].s << " && " << tgtBV.x.data[k].s << std::endl;
			//}
		}
	} else {
		//std::cerr << "pos=" << refBV.x.pos.size() << std::endl;
		for(int k = 0; k < tgtBV.x.pos.size(); ++k){
			const uint32_t reg = tgtBV.x.pos[k];
			if(!(refBV.x.data[reg].s && tgtBV.x.data[reg].s)) continue;
			//if(refBV.x.data[reg].s && tgtBV.x.data[reg].s){
				//uint32_t i = 0, j = 0;
				//uint32_t ii = 0, jj = 0;

				//std::cerr << "checking=" << tgtBV.x.data[reg].pos.size() << std::endl;
				const twk_igt_list::ilist_cont& ref = refBV.x.data[reg];
				const twk_igt_list::ilist_cont& tgt = tgtBV.x.data[reg];
				const std::vector<uint32_t>& posA = ref.pos;
				const std::vector<uint32_t>& posB = tgt.pos;

				if(posA.size() < posB.size()){
					for(int i = 0; i < posA.size(); ++i){
						const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)ref.data[posA[i]].vals);
						const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)tgt.data[posA[i]].vals);
						popcnt128(helper.alleleCounts[TWK_LD_ALTALT], _mm_and_si128(vectorA,vectorB));
					}
				} else {
					for(int i = 0; i < posB.size(); ++i){
						const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)ref.data[posB[i]].vals);
						const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)tgt.data[posB[i]].vals);
						popcnt128(helper.alleleCounts[TWK_LD_ALTALT], _mm_and_si128(vectorA,vectorB));
					}
				}

				/*while(true){
					// select
					const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)refBV.x.data[reg].data[i].vals);
					const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)tgtBV.x.data[reg].data[j].vals);

					popcnt128(helper.alleleCounts[TWK_LD_ALTALT],
							_mm_and_si128(
							  _mm_and_si128(vectorA,vectorB),
											_mm_set1_epi32(mask_lookup[refBV.x.data[reg].data[i].bin == tgtBV.x.data[reg].data[j].bin])));

					//std::cerr << "Diff@" << k << " -> " << refBV.x.data[k].data[i].bin << " and " << tgtBV.x.data[k].data[j].bin << std::endl;
					ii = i; jj = j;
					i += (refBV.x.data[reg].data[ii].bin == tgtBV.x.data[reg].data[jj].bin); // A == B
					j += (refBV.x.data[reg].data[ii].bin == tgtBV.x.data[reg].data[jj].bin); // A == B
					j += refBV.x.data[reg].data[ii].bin > tgtBV.x.data[reg].data[jj].bin; // A > B
					i += refBV.x.data[reg].data[ii].bin < tgtBV.x.data[reg].data[jj].bin; // A < B

					//std::cerr << i << "/" << ref.n << " and " << j << "/" << tgt.n << std::endl;
					if(i == refBV.x.data[reg].n || j == tgtBV.x.data[reg].n) break;
				}*\/
				//std::cerr << "overlap=" << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
			//}
			//else {
				//std::cerr << "no overlap = " << refBV.x.data[k].s << " && " << tgtBV.x.data[k].s << std::endl;
			//}
		}
	}

	/*uint32_t i = 0, j = 0;
	const uint32_t mask_lookup[2] = {(uint32_t)0, ~(uint32_t)0};
	uint32_t ii = 0, jj = 0;
	while(true){
		// select
		const VECTOR_TYPE vectorA = _mm_load_si128((const VECTOR_TYPE*)ref.data[i].vals);
		const VECTOR_TYPE vectorB = _mm_load_si128((const VECTOR_TYPE*)tgt.data[j].vals);

		popcnt128(helper.alleleCounts[TWK_LD_ALTALT],
				_mm_and_si128(
				  _mm_and_si128(vectorA,vectorB),
				                _mm_set1_epi32(mask_lookup[ref.data[i].bin == tgt.data[j].bin])));

		//std::cerr << "Diff=" << ref.data[i].bin << " and " << tgt.data[j].bin << std::endl;
		ii = i; jj = j;
		i += (ref.data[ii].bin == tgt.data[jj].bin); // A == B
		j += (ref.data[ii].bin == tgt.data[jj].bin); // A == B
		j += ref.data[ii].bin > tgt.data[jj].bin; // A > B
		i += ref.data[ii].bin < tgt.data[jj].bin; // A < B

		//std::cerr << i << "/" << ref.n << " and " << j << "/" << tgt.n << std::endl;
		if(i == ref.n || j == tgt.n) break;
	}*\/

#if TWK_SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	//const uint32_t n_list = refBV.l_list < tgtBV.l_list ? refBV.l_list : tgtBV.l_list;
	const uint32_t n_list = refBV.x_bins < tgtBV.x_bins ? refBV.x_bins : tgtBV.x_bins;
	perf->cycles[n_list] += ticks_per_iter.count();
	++perf->freq[n_list];
#endif

	helper.alleleCounts[TWK_LD_ALTREF] = b1.list[p1].l_list - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFALT] = b2.list[p2].l_list - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFREF] = 2*n_samples - ((b1.list[p1].l_list + b2.list[p2].l_list) - helper.alleleCounts[TWK_LD_ALTALT]);
	++n_method[0];

#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "listS=" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
#endif

#if TWK_SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}
*/

bool twk_ld_engine::PhasedBitmap(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
	helper.ResetPhased();
	const twk1_ldd_blk::bitmap_type& refB = b1.bitmap[p1];
	const twk1_ldd_blk::bitmap_type& tgtB = b2.bitmap[p2];
	// Debug timings
	#if TWK_SLAVE_DEBUG_MODE == 1
		typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
		auto t0 = std::chrono::high_resolution_clock::now();
	#endif
	helper.alleleCounts[TWK_LD_ALTALT] = refB.logicalandcount(tgtB);
	helper.alleleCounts[TWK_LD_ALTREF] = b1.blk->rcds[p1].ac - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFALT] = b2.blk->rcds[p2].ac - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFREF] = 2*n_samples - ((b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac) - helper.alleleCounts[TWK_LD_ALTALT]);
	++n_method[1];

#if TWK_SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	const uint32_t n_cycles = b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac;
	perf->cycles[n_cycles] += ticks_per_iter.count();
	++perf->freq[n_cycles];
#endif

#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "bitmap=" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
#endif

#if TWK_SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::PhasedVectorized(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
#if TWK_SLAVE_DEBUG_MODE == 0
	if(b1.blk->rcds[p1].gt_missing == false && b2.blk->rcds[p2].gt_missing == false){
		return(this->PhasedVectorizedNoMissing(b1,p1,b2,p2,perf));
	}
#endif

	helper.ResetPhased();
	helper_simd.counters[0] = 0;
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;
	// Data
	const twk_igt_vec& block1 = b1.vec[p1];
	const twk_igt_vec& block2 = b2.vec[p2];
	const uint64_t* const arrayA = (const uint64_t* const)block1.data;
	const uint64_t* const arrayB = (const uint64_t* const)block2.data;
	const uint64_t* const arrayA_mask = b1.blk->rcds[p1].gt_missing ? (const uint64_t* const)block1.mask : (const uint64_t* const)mask_placeholder;
	const uint64_t* const arrayB_mask = b2.blk->rcds[p2].gt_missing ? (const uint64_t* const)block2.mask : (const uint64_t* const)mask_placeholder;

#if SIMD_AVAILABLE == 1
	const uint32_t frontSmallest = block1.front_zero < block2.front_zero ? block1.front_zero : block2.front_zero;
	const uint32_t tailSmallest  = block1.tail_zero  < block2.tail_zero  ? block1.tail_zero  : block2.tail_zero;
	const uint32_t frontBonus    = block1.front_zero != frontSmallest ? block1.front_zero : block2.front_zero;
	const uint32_t tailBonus     = block1.tail_zero  != tailSmallest  ? block1.tail_zero  : block2.tail_zero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	const VECTOR_TYPE* const vectorA_mask = (const VECTOR_TYPE* const)arrayA_mask;
	const VECTOR_TYPE* const vectorB_mask = (const VECTOR_TYPE* const)arrayB_mask;
	VECTOR_TYPE a, b;

	VECTOR_TYPE __intermediate, masks;

	uint32_t i = frontSmallest;

// Debug timings
#if TWK_SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_SHORT {                                                     \
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);              \
	a = _mm_load_si128( (__m128i*) &vectorA[i] ); \
	b = _mm_load_si128( (__m128i*) &vectorB[i] ); \
	__intermediate  = PHASED_REFREF_MASK(a, b, masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_REFREF], __intermediate);       \
	__intermediate  = PHASED_ALTREF_MASK(a, b, masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_ALTREF], __intermediate);       \
	__intermediate  = PHASED_REFALT_MASK(a, b, masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_REFALT], __intermediate);       \
	i += 1;                                                              \
}

#define ITER {                                                           \
	masks   = MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);              \
	__m128i a = _mm_load_si128( (__m128i*) &vectorA[i] ); \
	__m128i b = _mm_load_si128( (__m128i*) &vectorB[i] ); \
	__intermediate  = PHASED_ALTALT_MASK(a, b, masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_ALTALT], __intermediate);       \
	__intermediate  = PHASED_REFREF_MASK(a, b, masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_REFREF], __intermediate);       \
	__intermediate  = PHASED_ALTREF_MASK(a, b, masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_ALTREF], __intermediate);       \
	__intermediate  = PHASED_REFALT_MASK(a, b, masks); \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_REFALT], __intermediate);       \
	i += 1;                                                              \
}

	for( ; i < frontBonus; ) 					  	 ITER_SHORT // Not possible to be ALT-ALT
	for( ; i < this->vector_cycles - tailBonus; )  	 ITER
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT // Not possible to be ALT-ALT

#undef ITER
#undef ITER_SHORT
	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif

	uint64_t b_altalt, b_refref, b_refalt, b_altref;
	for(; k < this->byte_width; ++k){
		b_altalt  = (arrayA[k] & arrayB[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_refref  = ((~arrayA[k]) & (~arrayB[k])) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_altref  = ((arrayA[k] ^ arrayB[k]) & arrayA[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_refalt  = ((arrayA[k] ^ arrayB[k]) & arrayB[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);
		helper_simd.counters[TWK_LD_SIMD_REFREF] += POPCOUNT_ITER(b_refalt);
		helper_simd.counters[TWK_LD_SIMD_REFALT] += POPCOUNT_ITER(b_refalt);
		helper_simd.counters[TWK_LD_SIMD_ALTREF] += POPCOUNT_ITER(b_altref);
		helper_simd.counters[TWK_LD_SIMD_ALTALT] += POPCOUNT_ITER(b_altalt);
	}

	helper.alleleCounts[TWK_LD_REFALT] = helper_simd.counters[TWK_LD_SIMD_ALTREF];
	helper.alleleCounts[TWK_LD_ALTREF] = helper_simd.counters[TWK_LD_SIMD_REFALT];
	helper.alleleCounts[TWK_LD_ALTALT] = helper_simd.counters[TWK_LD_SIMD_ALTALT];
	helper.alleleCounts[TWK_LD_REFREF] = helper_simd.counters[TWK_LD_SIMD_REFREF] + (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT*2 - phased_unbalanced_adjustment;
	//helper.haplotypeCounts[0] += (tailSmallest + frontSmallest) * GENOTYPE_TRIP_COUNT;
	++n_method[2];

#if TWK_SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
	//std::cout << "V\t" << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif


	//this->setFLAGs(block1, block2);
	//return(this->CalculateLDPhasedMath());
	//std::cerr << "m1front=" << frontBonus << "," << tailBonus << "," << phased_unbalanced_adjustment << std::endl;
#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "m1 " << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
#endif

#if TWK_SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::PhasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
	helper.ResetPhased();
	helper_simd.counters[TWK_LD_SIMD_ALTALT] = 0;

	const twk_igt_vec& block1 = b1.vec[p1];
	const twk_igt_vec& block2 = b2.vec[p2];
	const uint64_t* const arrayA = (const uint64_t* const)block1.data;
	const uint64_t* const arrayB = (const uint64_t* const)block2.data;

	// Debug timings
#if TWK_SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#if SIMD_AVAILABLE == 1
	const uint32_t frontSmallest = block1.front_zero < block2.front_zero ? block1.front_zero : block2.front_zero;
	const uint32_t tailSmallest  = block1.tail_zero  < block2.tail_zero  ? block1.tail_zero  : block2.tail_zero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	VECTOR_TYPE __intermediate, a, b;

#define ITER_SHORT {                                              \
	a = _mm_load_si128( (__m128i*) &vectorA[i] ); \
	b = _mm_load_si128( (__m128i*) &vectorB[i] ); \
	__intermediate  = PHASED_ALTALT(a, b);      \
	popcnt128(helper_simd.counters[TWK_LD_SIMD_ALTALT], __intermediate);\
	i += 1;                                                       \
}

	uint32_t i = frontSmallest;
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT


#undef ITER_SHORT
	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif

	uint64_t b_altalt, b_refref, b_refalt, b_altref;
	for(; k < this->byte_width; ++k){
		b_altalt  = (arrayA[k] & arrayB[k]);
		helper_simd.counters[TWK_LD_SIMD_ALTALT] += POPCOUNT_ITER(b_altalt);
	}
	helper.alleleCounts[TWK_LD_ALTALT] = helper_simd.counters[TWK_LD_SIMD_ALTALT];
	helper.alleleCounts[TWK_LD_ALTREF] = b1.blk->rcds[p1].ac - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFALT] = b2.blk->rcds[p2].ac - helper.alleleCounts[TWK_LD_ALTALT];
	helper.alleleCounts[TWK_LD_REFREF] = 2*n_samples - ((b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac) - helper.alleleCounts[TWK_LD_ALTALT]);
	++n_method[3];

#if TWK_SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
	//std::cout << "V\t" << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << "\t" << ticks_per_iter.count() << '\n';
#endif

#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "m3=" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
#endif

	//this->setFLAGs(block1, block2);
	//return(this->CalculateLDPhasedMathSimple(block1, block2));
#if TWK_SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::UnphasedVectorized(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
#if TWK_SLAVE_DEBUG_MODE == 0
	if(b1.blk->rcds[p1].gt_missing == false && b2.blk->rcds[p2].gt_missing == false){
		return(this->UnphasedVectorizedNoMissing(b1,p1,b2,p2,perf));
	}
#endif

	helper.ResetUnphased();

	helper_simd.counters[0] = 0;
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;
	helper_simd.counters[4] = 0;
	helper_simd.counters[5] = 0;
	helper_simd.counters[6] = 0;
	helper_simd.counters[7] = 0;
	helper_simd.counters[8] = 0;
	helper_simd.counters[9] = 0;

	// Data
	const twk_igt_vec& block1 = b1.vec[p1];
	const twk_igt_vec& block2 = b2.vec[p2];
	const uint64_t* const arrayA = (const uint64_t* const)block1.data;
	const uint64_t* const arrayB = (const uint64_t* const)block2.data;
	const uint64_t* const arrayA_mask = b1.blk->rcds[p1].gt_missing ? (const uint64_t* const)block1.mask : (const uint64_t* const)mask_placeholder;
	const uint64_t* const arrayB_mask = b2.blk->rcds[p2].gt_missing ? (const uint64_t* const)block2.mask : (const uint64_t* const)mask_placeholder;

#if SIMD_AVAILABLE == 1
	const uint32_t frontSmallest = block1.front_zero < block2.front_zero ? block1.front_zero : block2.front_zero;
	const uint32_t tailSmallest  = block1.tail_zero  < block2.tail_zero  ? block1.tail_zero  : block2.tail_zero;
	const uint32_t frontBonus    = block1.front_zero != frontSmallest ? block1.front_zero : block2.front_zero;
	const uint32_t tailBonus     = block1.tail_zero  != tailSmallest  ? block1.tail_zero  : block2.tail_zero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;
	const VECTOR_TYPE* const vectorA_mask = (const VECTOR_TYPE* const)arrayA_mask;
	const VECTOR_TYPE* const vectorB_mask = (const VECTOR_TYPE* const)arrayB_mask;

	VECTOR_TYPE __intermediate, mask, a, b;

	uint32_t i = frontSmallest;
	VECTOR_TYPE altalt, refref, altref, refalt;
	//VECTOR_TYPE t00, t05, t50, t55, combTop, combLeft, combRight, combBottom;


// Debug timings
#if TWK_SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_BASE {                                             \
	mask	= MASK_MERGE(vectorA_mask[i], vectorB_mask[i]);     \
	a = _mm_load_si128( (__m128i*) &vectorA[i] ); \
	b = _mm_load_si128( (__m128i*) &vectorB[i] ); \
	refref  = PHASED_REFREF_MASK(a, b, mask);	\
	refalt  = PHASED_REFALT_MASK(a, b, mask);	\
	altref  = PHASED_ALTREF_MASK(a, b, mask);	\
	altalt  = PHASED_ALTALT_MASK(a, b, mask);	\
	i += 1;                                                     \
}

#define ITER_SHORT {                                                        \
	ITER_BASE                                                               \
	__intermediate = FILTER_UNPHASED_SPECIAL(refref);                       \
	popcnt128(helper_simd.counters[0], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);    \
	popcnt128(helper_simd.counters[1], __intermediate);                      \
	__intermediate = FILTER_UNPHASED(altref, altref);                       \
	popcnt128(helper_simd.counters[2], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);    \
	popcnt128(helper_simd.counters[3], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refref, altalt, altalt, refref);  \
	popcnt128(helper_simd.counters[4], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altref, altref, refalt);  \
	popcnt128(helper_simd.counters[4], __intermediate);                      \
	__intermediate = FILTER_UNPHASED(refalt,refalt);                        \
	popcnt128(helper_simd.counters[6], __intermediate);                      \
}

#define ITER_LONG {                                                         \
	ITER_SHORT                                                              \
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref);    \
	popcnt128(helper_simd.counters[5], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt);  \
	popcnt128(helper_simd.counters[7], __intermediate);                      \
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt);                       \
	popcnt128(helper_simd.counters[8], __intermediate);                      \
}

	for( ; i < frontBonus; ) 					  	 ITER_SHORT
	for( ; i < this->vector_cycles - tailBonus; )  	 ITER_LONG
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT
	//for( ; i < this->vector_cycles; )  	 ITER_LONG

#undef ITER_LONG
#undef ITER_SHORT
#undef ITER_BASE

	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif

	//k=0;
	uint64_t b_altalt, b_refref, b_refalt, b_altref;
	for(; k < this->byte_width; ++k){
		b_altalt  = (arrayA[k] & arrayB[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_refref  = ((~arrayA[k]) & (~arrayB[k])) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_altref  = ((arrayA[k] ^ arrayB[k]) & arrayA[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);
		b_refalt  = ((arrayA[k] ^ arrayB[k]) & arrayB[k]) & ~(arrayA_mask[k] | arrayB_mask[k]);

		helper_simd.counters[0] += POPCOUNT_ITER(FILTER_UNPHASED_64_SPECIAL(b_refref));
		helper_simd.counters[1] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_refref, b_refalt, b_refalt, b_refref));
		helper_simd.counters[2] += POPCOUNT_ITER(FILTER_UNPHASED_64(b_refalt, b_refalt));
		helper_simd.counters[3] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_refref, b_altref, b_altref, b_refref));
		helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_refref, b_altalt, b_altalt, b_refref));
		helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_refalt, b_altref, b_altref, b_refalt));
		helper_simd.counters[5] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_refalt, b_altalt, b_altalt, b_refalt));
		helper_simd.counters[6] += POPCOUNT_ITER(FILTER_UNPHASED_64(b_altref, b_altref));
		helper_simd.counters[7] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_altref, b_altalt, b_altalt, b_altref));
		helper_simd.counters[8] += POPCOUNT_ITER(FILTER_UNPHASED_64_SPECIAL(b_altalt));
	}

	//std::cerr << "0 count=" << helper_simd.counters[0] << " adjust=" << unphased_unbalanced_adjustment << std::endl;
	helper.alleleCounts[TWK_LD_REFREF]  = helper_simd.counters[0] - this->unphased_unbalanced_adjustment;
	helper.alleleCounts[TWK_LD_REFREF] += (frontSmallest + tailSmallest) * GENOTYPE_TRIP_COUNT;
	helper.alleleCounts[TWK_LD_REFALT]  = helper_simd.counters[1];
	helper.alleleCounts[TWK_LD_ALTALT]  = helper_simd.counters[2];
	helper.alleleCounts[16] = helper_simd.counters[3];
	helper.alleleCounts[17] = helper_simd.counters[4];
	helper.alleleCounts[21] = helper_simd.counters[5];
	helper.alleleCounts[80] = helper_simd.counters[6];
	helper.alleleCounts[81] = helper_simd.counters[7];
	helper.alleleCounts[85] = helper_simd.counters[8];
	++n_method[4];

#if TWK_SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	//std::cout << ticks_per_iter.count() << '\n';
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
#endif

#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "vum =" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTALT]
	          << "," << helper.alleleCounts[16] << "," << helper.alleleCounts[17] << "," << helper.alleleCounts[21]
	          << "," << helper.alleleCounts[80] << "," << helper.alleleCounts[81] << "," << helper.alleleCounts[85] << std::endl;
#endif

#if TWK_SLAVE_DEBUG_MODE != 1
	return(UnphasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::UnphasedVectorizedNoMissing(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
	helper.ResetUnphased();

	helper_simd.counters[0] = 0;
	helper_simd.counters[1] = 0;
	helper_simd.counters[2] = 0;
	helper_simd.counters[3] = 0;
	helper_simd.counters[4] = 0;
	helper_simd.counters[5] = 0;
	helper_simd.counters[6] = 0;
	helper_simd.counters[7] = 0;
	helper_simd.counters[8] = 0;
	helper_simd.counters[9] = 0;

	// Data
	const twk_igt_vec& block1 = b1.vec[p1];
	const twk_igt_vec& block2 = b2.vec[p2];
	const uint64_t* const arrayA = (const uint64_t* const)block1.data;
	const uint64_t* const arrayB = (const uint64_t* const)block2.data;

#if SIMD_AVAILABLE == 1
	const uint32_t frontSmallest = block1.front_zero < block2.front_zero ? block1.front_zero : block2.front_zero;
	int i = frontSmallest;
	const uint32_t tailSmallest  = block1.tail_zero  < block2.tail_zero  ? block1.tail_zero  : block2.tail_zero;
	const uint32_t frontBonus    = block1.front_zero != frontSmallest ? block1.front_zero : block2.front_zero;
	const uint32_t tailBonus     = block1.tail_zero  != tailSmallest  ? block1.tail_zero  : block2.tail_zero;

	const VECTOR_TYPE* const vectorA = (const VECTOR_TYPE* const)arrayA;
	const VECTOR_TYPE* const vectorB = (const VECTOR_TYPE* const)arrayB;

	VECTOR_TYPE altalt, refref, altref, refalt;
	VECTOR_TYPE __intermediate, a, b;

// Debug timings
#if TWK_SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

#define ITER_BASE {                               \
	a = _mm_load_si128( (__m128i*) &vectorA[i] ); \
	b = _mm_load_si128( (__m128i*) &vectorB[i] ); \
	refref  = PHASED_REFREF(a, b);                \
	refalt  = PHASED_REFALT(a, b);                \
	altref  = PHASED_ALTREF(a, b);                \
	altalt  = PHASED_ALTALT(a, b);                \
	i += 1;                                       \
}

#define ITER_SHORT {                                                      \
	ITER_BASE                                                             \
	__intermediate = FILTER_UNPHASED_SPECIAL(refref);                     \
	popcnt128(helper_simd.counters[0], __intermediate);                   \
	__intermediate = FILTER_UNPHASED_PAIR(refref,altref, altref,refref);  \
	popcnt128(helper_simd.counters[1], __intermediate);                   \
	__intermediate = FILTER_UNPHASED(altref, altref);                     \
	popcnt128(helper_simd.counters[2], __intermediate);                   \
	__intermediate = FILTER_UNPHASED_PAIR(refref,refalt, refalt,refref);  \
	popcnt128(helper_simd.counters[3], __intermediate);                   \
	__intermediate = FILTER_UNPHASED(refalt,refalt);                      \
	popcnt128(helper_simd.counters[6], __intermediate);                   \
}

#define ITER_LONG {                                                        \
	ITER_SHORT                                                             \
	__intermediate = FILTER_UNPHASED_PAIR(altref,altalt, altalt,altref);   \
	popcnt128(helper_simd.counters[5], __intermediate);                    \
	__intermediate = FILTER_UNPHASED_PAIR(refalt, altalt, altalt, refalt); \
	popcnt128(helper_simd.counters[7], __intermediate);                    \
	__intermediate = FILTER_UNPHASED_SPECIAL(altalt);                      \
	popcnt128(helper_simd.counters[8], __intermediate);                    \
}

	for( ; i < frontBonus; )                         ITER_SHORT
	for( ; i < this->vector_cycles - tailBonus; )    ITER_LONG
	for( ; i < this->vector_cycles - tailSmallest; ) ITER_SHORT

#undef ITER_LONG
#undef ITER_SHORT
#undef ITER_BASE

	uint32_t k = this->byte_aligned_end;
#else
	uint32_t k = 0;
#endif

	//k = 0;
	//std::cerr << "at=" << byte_aligned_end << "/" << byte_width << " " << "vector cycles=" << vector_cycles << std::endl;
	uint64_t b_altalt, b_refref, b_refalt, b_altref;
	for(; k < this->byte_width; ++k){
		b_altalt  = (arrayA[k] & arrayB[k]);
		b_refref  = ((~arrayA[k]) & (~arrayB[k]));
		b_altref  = ((arrayA[k] ^ arrayB[k]) & arrayA[k]);
		b_refalt  = ((arrayA[k] ^ arrayB[k]) & arrayB[k]);

		helper_simd.counters[0] += POPCOUNT_ITER(FILTER_UNPHASED_64_SPECIAL(b_refref));
		helper_simd.counters[1] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_refref, b_refalt, b_refalt, b_refref));
		helper_simd.counters[2] += POPCOUNT_ITER(FILTER_UNPHASED_64(b_refalt, b_refalt));
		helper_simd.counters[3] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_refref, b_altref, b_altref, b_refref));
		//helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refref, b_altalt, b_altalt, b_refref));
		//helper_simd.counters[4] += POPCOUNT_ITER(FILTER_UNPHASED_BYTE_PAIR(b_refalt, b_altref, b_altref, b_refalt));

		helper_simd.counters[5] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_refalt, b_altalt, b_altalt, b_refalt));
		helper_simd.counters[6] += POPCOUNT_ITER(FILTER_UNPHASED_64(b_altref, b_altref));
		helper_simd.counters[7] += POPCOUNT_ITER(FILTER_UNPHASED_64_PAIR(b_altref, b_altalt, b_altalt, b_altref));
		helper_simd.counters[8] += POPCOUNT_ITER(FILTER_UNPHASED_64_SPECIAL(b_altalt));
	}

	helper.alleleCounts[TWK_LD_REFREF]  = helper_simd.counters[0] - this->unphased_unbalanced_adjustment;
	helper.alleleCounts[TWK_LD_REFREF] += (frontSmallest + tailSmallest) * GENOTYPE_TRIP_COUNT;
	helper.alleleCounts[TWK_LD_REFALT]  = helper_simd.counters[1];
	helper.alleleCounts[TWK_LD_ALTALT]  = helper_simd.counters[2];
	helper.alleleCounts[16] = helper_simd.counters[3];
	helper.alleleCounts[21] = helper_simd.counters[5];
	helper.alleleCounts[80] = helper_simd.counters[6];
	helper.alleleCounts[81] = helper_simd.counters[7];
	helper.alleleCounts[85] = helper_simd.counters[8];
	//helper[17] = helper_simd.counters[4];
	helper.alleleCounts[17] = n_samples - (helper.alleleCounts[TWK_LD_REFREF] +  helper.alleleCounts[TWK_LD_REFALT] + helper.alleleCounts[TWK_LD_ALTALT] + helper.alleleCounts[16] + helper.alleleCounts[21] + helper.alleleCounts[80] + helper.alleleCounts[81] + helper.alleleCounts[85]);
	++n_method[5];

#if TWK_SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	//std::cout << ticks_per_iter.count() << '\n';
	perf->cycles[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac] += ticks_per_iter.count();
	++perf->freq[b1.blk->rcds[p1].ac + b2.blk->rcds[p2].ac];
#endif

#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "vu  =" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTALT]
		      << "," << helper.alleleCounts[16] << "," << helper.alleleCounts[17] << "," << helper.alleleCounts[21]
		      << "," << helper.alleleCounts[80] << "," << helper.alleleCounts[81] << "," << helper.alleleCounts[85] << std::endl;
#endif

	//this->setFLAGs(block1, block2);
#if TWK_SLAVE_DEBUG_MODE != 1
	return(UnphasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::PhasedRunlength(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
	helper.ResetPhased();
#if TWK_SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	const twk1_gt_t& gt1 = *b1.blk->rcds[p1].gt;
	const twk1_gt_t& gt2 = *b2.blk->rcds[p2].gt;

	uint32_t lenA = gt1.GetLength(0);
	uint32_t lenB = gt2.GetLength(0);
	uint8_t currentMixL = (gt1.GetRefA(0) << 2) | gt2.GetRefA(0);
	uint8_t currentMixR = (gt1.GetRefB(0) << 2) | gt2.GetRefB(0);
	uint32_t offsetA = 0;
	uint32_t offsetB = 0;
	uint32_t add;
	// Branchless
	//__m128* voffset = reinterpret_cast<__m128*>(&offsetA);
	//__m128* vmix    = reinterpret_cast<__m128*>(&currentMixL);
	//__m128* vlen    = reinterpret_cast<__m128*>(&lenA);

	//_mm_blendv_ps(voffest)

	// tests
	// A = lenA > lenB ->  _mm_cmpgt_epi32
	// B = lenB < lenA ->  _mm_cmplt_epi32
	// else if merge above -> _mm_and_ps, _mm_cmpeq_epi32
	// A & B == 0

	while(true){
		if(lenA > lenB){ // If processed run length A > processed run length B
			lenA -= lenB;
			add = lenB;
			lenB = gt2.GetLength(++offsetB);
		} else if(lenA < lenB){ // If processed run length A < processed run length B
			lenB -= lenA;
			add = lenA;
			lenA = gt1.GetLength(++offsetA);
		} else { // If processed run length A == processed run length B
			add = lenB;
			lenA = gt1.GetLength(++offsetA);
			lenB = gt2.GetLength(++offsetB);
		}
		helper.alleleCounts[currentMixL] += add;
		helper.alleleCounts[currentMixR] += add;

		// Exit condition
		if(offsetA == gt1.n || offsetB == gt2.n){
			if(offsetA != gt1.n || offsetB != gt2.n){
				std::cerr << utility::timestamp("FATAL") << "Failed to exit equally!\n" << offsetA << "/" << gt1.n << " and " << offsetB << "/" << gt2.n << std::endl;
				exit(1);
			}
			break;
		}

		currentMixL = (gt1.GetRefA(offsetA) << 2) | gt2.GetRefA(offsetB);
		currentMixR = (gt1.GetRefB(offsetA) << 2) | gt2.GetRefB(offsetB);
	}
	++n_method[6];

#if TWK_SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[gt1.n + gt2.n] += ticks_per_iter.count();
	++perf->freq[gt1.n + gt2.n];
	//std::cout << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif

#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "rleP=" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
#endif

	//std::cerr << "p= " << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT] << "," << helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT] << std::endl;

#if TWK_SLAVE_DEBUG_MODE != 1
	return(PhasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::UnphasedRunlength(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2, twk_ld_perf* perf){
	helper.ResetUnphased();
#if TWK_SLAVE_DEBUG_MODE == 1
	typedef std::chrono::duration<double, typename std::chrono::high_resolution_clock::period> Cycle;
	auto t0 = std::chrono::high_resolution_clock::now();
#endif

	const twk1_gt_t& gt1 = *b1.blk->rcds[p1].gt;
	const twk1_gt_t& gt2 = *b2.blk->rcds[p2].gt;

	uint32_t lenA = gt1.GetLength(0);
	uint32_t lenB = gt2.GetLength(0);
	uint8_t  currentMix = (gt1.GetRefA(0) << 6) | (gt1.GetRefB(0) << 4) | (gt2.GetRefA(0) << 2) | (gt2.GetRefB(0));
	uint32_t offsetA = 0;
	uint32_t offsetB = 0;
	uint32_t add;

	while(true){
		if(lenA > lenB){ // If processed run length A > processed run length B
			lenA -= lenB;
			add = lenB;
			lenB = gt2.GetLength(++offsetB);
		} else if(lenA < lenB){ // If processed run length A < processed run length B
			lenB -= lenA;
			add = lenA;
			lenA = gt1.GetLength(++offsetA);
		} else { // If processed run length A == processed run length B
			add = lenB;
			lenA = gt1.GetLength(++offsetA);
			lenB = gt2.GetLength(++offsetB);
		}
		helper.alleleCounts[currentMix] += add;
		//std::cerr << "adding: " << std::bitset<8>(currentMix) << std::endl;
		assert(currentMix < 171);

		// Exit condition
		if(offsetA == gt1.n && offsetB == gt2.n){
			if(offsetA != gt1.n || offsetB != gt2.n){
				std::cerr << utility::timestamp("FATAL") << "Failed to exit equally!\n" << offsetA << "/" << gt1.n << " and " << offsetB << "/" << gt2.n << std::endl;
				exit(1);
			}
			break;
		}

		currentMix = (gt1.GetRefA(offsetA) << 6) | (gt1.GetRefB(offsetA) << 4) | (gt2.GetRefA(offsetB) << 2) | (gt2.GetRefB(offsetB));
	}
	++n_method[7];

#if TWK_SLAVE_DEBUG_MODE == 1
	auto t1 = std::chrono::high_resolution_clock::now();
	auto ticks_per_iter = Cycle(t1-t0);
	perf->cycles[gt1.n + gt2.n] += ticks_per_iter.count();
	++perf->freq[gt1.n + gt2.n];
	//std::cout << a.getMeta().MAF*this->samples + b.getMeta().MAF*this->samples << '\t' << ticks_per_iter.count() << '\n';
#endif

#if TWK_SLAVE_DEBUG_MODE == 2
	std::cerr << "rleU=" << helper.alleleCounts[TWK_LD_REFREF] << "," << helper.alleleCounts[TWK_LD_REFALT]+ helper.alleleCounts[TWK_LD_ALTREF] << "," << helper.alleleCounts[TWK_LD_ALTALT]
						 << "," << helper.alleleCounts[16]+helper.alleleCounts[64] << "," << helper.alleleCounts[17]+ helper.alleleCounts[20]+ helper.alleleCounts[65] + helper.alleleCounts[68] << "," << helper.alleleCounts[21]+ helper.alleleCounts[69]
						 << "," << helper.alleleCounts[80] << "," << helper.alleleCounts[81]+ helper.alleleCounts[84] << "," << helper.alleleCounts[85] << std::endl;
#endif

#if TWK_SLAVE_DEBUG_MODE != 1
	return(UnphasedMath(b1,p1,b2,p2));
#else
	return(true);
#endif
}

bool twk_ld_engine::PhasedMath(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2){
	// Total amount of non-missing alleles
	helper.totalHaplotypeCounts = helper.alleleCounts[TWK_LD_REFREF] + helper.alleleCounts[TWK_LD_REFALT] +
								  helper.alleleCounts[TWK_LD_ALTREF] + helper.alleleCounts[TWK_LD_ALTALT];

	// All values are missing
	if(helper.totalHaplotypeCounts < TWK_MINIMUM_ALLOWED_ALLELES){
		//++this->insufficent_alleles;
		return false;
	}

	// Todo: temporary
	if(helper.alleleCounts[TWK_LD_REFREF] < helper.alleleCounts[TWK_LD_ALTALT]){
		if(helper.alleleCounts[TWK_LD_REFALT] + helper.alleleCounts[TWK_LD_ALTREF] + helper.alleleCounts[TWK_LD_REFREF] < 5){
			//std::cerr << "return=" << helper.alleleCounts[TWK_LD_REFALT] + helper.alleleCounts[TWK_LD_ALTREF] + helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
			cur_rcd.controller = 0;
			return false;
		}
	} else {
		if(helper.alleleCounts[TWK_LD_ALTALT] + helper.alleleCounts[TWK_LD_REFALT] + helper.alleleCounts[TWK_LD_ALTREF] < 5){
			//std::cerr << "return=" << helper.alleleCounts[TWK_LD_REFALT] + helper.alleleCounts[TWK_LD_ALTREF] + helper.alleleCounts[TWK_LD_ALTALT] << std::endl;
			cur_rcd.controller = 0;
			return false;
		}
	}

	// Haplotype frequencies
	double pA = TWK_HAP_FREQ(helper,0); // pA
	double qA = TWK_HAP_FREQ(helper,1); // qA
	double pB = TWK_HAP_FREQ(helper,4); // pB
	double qB = TWK_HAP_FREQ(helper,5); // qB

	if(pA*qB - qA*pB == 0){ return false; }

	// Allelic frequencies
	const double g0 = ((double)helper.alleleCounts[TWK_LD_REFREF] + helper.alleleCounts[TWK_LD_REFALT]) / (helper.totalHaplotypeCounts);
	const double g1 = ((double)helper.alleleCounts[TWK_LD_ALTREF] + helper.alleleCounts[TWK_LD_ALTALT]) / (helper.totalHaplotypeCounts);
	const double h0 = ((double)helper.alleleCounts[TWK_LD_REFREF] + helper.alleleCounts[TWK_LD_ALTREF]) / (helper.totalHaplotypeCounts);
	const double h1 = ((double)helper.alleleCounts[TWK_LD_REFALT] + helper.alleleCounts[TWK_LD_ALTALT]) / (helper.totalHaplotypeCounts);

	cur_rcd.D = pA*qB - qA*pB;
	cur_rcd.R2 = cur_rcd.D*cur_rcd.D / (g0*g1*h0*h1);
	if(cur_rcd.R2 < settings.minR2 || cur_rcd.R2 > settings.maxR2){
		cur_rcd.controller = 0;
		return false;
	}

	double dmax = 0;
	if(cur_rcd.D >= 0) dmax = g0*h1 < h0*g1 ? g0*h1 : h0*g1;
	else dmax = g0*g1 < h0*h1 ? -g0*g1 : -h0*h1;

	cur_rcd.Dprime = cur_rcd.D / dmax;

	if(cur_rcd.Dprime < settings.minDprime || cur_rcd.Dprime > settings.maxDprime){
		cur_rcd.controller = 0;
		return false;
	}

	// Calculate Fisher's exact test P-value
	double left,right,both;
	kt_fisher_exact(helper.alleleCounts[TWK_LD_REFREF],
					helper.alleleCounts[TWK_LD_REFALT],
					helper.alleleCounts[TWK_LD_ALTREF],
					helper.alleleCounts[TWK_LD_ALTALT],
					&left, &right, &both);

	if(both > settings.minP){
		cur_rcd.controller = 0;
		return false;
	}
	cur_rcd.P = both;
	cur_rcd.R = sqrt(cur_rcd.R2);

	cur_rcd.Apos = b1.blk->rcds[p1].pos;
	cur_rcd.Bpos = b2.blk->rcds[p2].pos;
	cur_rcd.ridA = b1.blk->rcds[p1].rid;
	cur_rcd.ridB = b2.blk->rcds[p2].rid;
	cur_rcd[TWK_LD_SIMD_REFREF] = helper.alleleCounts[TWK_LD_REFREF];
	cur_rcd[TWK_LD_SIMD_REFALT] = helper.alleleCounts[TWK_LD_REFALT];
	cur_rcd[TWK_LD_SIMD_ALTREF] = helper.alleleCounts[TWK_LD_ALTREF];
	cur_rcd[TWK_LD_SIMD_ALTALT] = helper.alleleCounts[TWK_LD_ALTALT];

	cur_rcd.SetLowACA(b1.blk->rcds[p1].ac < TWK_LOW_AC_THRESHOLD);
	cur_rcd.SetLowACB(b2.blk->rcds[p2].ac < TWK_LOW_AC_THRESHOLD);
	cur_rcd.SetCompleteLD(helper.alleleCounts[TWK_LD_REFREF] < 1 || helper.alleleCounts[TWK_LD_REFALT] < 1 || helper.alleleCounts[TWK_LD_ALTREF] < 1 || helper.alleleCounts[TWK_LD_ALTALT] < 1);
	cur_rcd.SetPerfectLD(cur_rcd.R2 > 0.99);
	cur_rcd.SetHasMissingValuesA(b1.blk->rcds[p1].an);
	cur_rcd.SetHasMissingValuesB(b2.blk->rcds[p2].an);
	int32_t diff = (int32_t)b1.blk->rcds[p1].pos - b2.blk->rcds[p2].pos;
	cur_rcd.SetLongRange(abs(diff) > TWK_LONG_RANGE_THRESHOLD && (cur_rcd.ridA == cur_rcd.ridB));
	cur_rcd.SetUsedPhasedMath();
	cur_rcd.SetSameContig(cur_rcd.ridA == cur_rcd.ridB);
	cur_rcd.SetInvalidHWEA(b1.blk->rcds[p1].hwe < TWK_INVALID_HWE_THRESHOLD);
	cur_rcd.SetInvalidHWEB(b2.blk->rcds[p2].hwe < TWK_INVALID_HWE_THRESHOLD);

	// Calculate Chi-Sq CV from 2x2 contingency table
	cur_rcd.ChiSqModel = 0;
	cur_rcd.ChiSqFisher = helper.totalHaplotypeCounts * cur_rcd.R2;

#if TWK_SLAVE_DEBUG_MODE == 3
	if(cur_rcd.R2 > 0.1){
		twk_debug_pos1 = b1.blk->rcds[p1].pos;
		twk_debug_pos1_2 = b2.blk->rcds[p2].pos;
//std::cout << "P\t" << b1.blk->rcds[p1].pos << "\t" << b2.blk->rcds[p2].pos << "\t" << helper.D << "\t" << helper.Dprime << "\t" << helper.R << "\t" << helper.R2 << '\n';
	}
#endif
	// If the number of rcds written is equal to the flush limit then
	// compress and write output.
	if(blk_f.n == n_lim || irecF.rid != cur_rcd.ridA || irecR.rid != cur_rcd.ridB){
		if(this->CompressBlock() == false) return false;
		irecF.rid    = cur_rcd.ridA;
		irecF.ridB   = cur_rcd.ridB;
		irecF.minpos = cur_rcd.Apos;
		irecF.maxpos = cur_rcd.Apos;
		irecR.rid    = cur_rcd.ridB;
		irecR.ridB   = cur_rcd.ridA;
		irecR.minpos = cur_rcd.Bpos;
		irecR.maxpos = cur_rcd.Bpos;
		//n_out_tick = 0;
	}
	// If ridB is mixed then set index value to -1.
	if(irecF.ridB != cur_rcd.ridB) irecF.ridB = -1;
	if(irecR.ridB != cur_rcd.ridA) irecR.ridB = -1;

	// Update index
	irecF.maxpos = cur_rcd.Apos;
	irecR.maxpos = cur_rcd.Bpos;

	// Add forward.
	blk_f += cur_rcd;
	// Swap tuple (ridA:posA) with (ridB:posB).
	uint32_t temp = cur_rcd.Apos;
	cur_rcd.Apos = cur_rcd.Bpos;
	cur_rcd.Bpos = temp;
	std::swap(cur_rcd.ridA,cur_rcd.ridB);
	// Add reverse.
	blk_r += cur_rcd;
	// Update tickers.
	n_out += 2; //n_out_tick += 2;

	// Update ticker.
	//if(n_out_tick == 300){ progress->n_out += 300; n_out_tick = 0; }
	//progress->n_out += 2;

	// Reset controller.
	cur_rcd.controller = 0;

	return true;
}

bool twk_ld_engine::UnphasedMath(const twk1_ldd_blk& b1, const uint32_t p1, const twk1_ldd_blk& b2, const uint32_t p2){
	// Total amount of non-missing alleles
	helper.totalHaplotypeCounts =
		  helper.alleleCounts[0]  + helper.alleleCounts[1]  + helper.alleleCounts[4]  + helper.alleleCounts[5]
	    + helper.alleleCounts[16] + helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[21]
	    + helper.alleleCounts[64] + helper.alleleCounts[65] + helper.alleleCounts[68] + helper.alleleCounts[69]
	    + helper.alleleCounts[80] + helper.alleleCounts[81] + helper.alleleCounts[84] + helper.alleleCounts[85];

	// All values are missing or too few
	if(helper.totalHaplotypeCounts < TWK_MINIMUM_ALLOWED_ALLELES){
		//++this->insufficent_alleles;
		return false;
	}

	// How many hets-hets is there? 0/1 0/1 or 1/0 1/0 or equivalent
	const uint64_t number_of_hets = helper.alleleCounts[17] + helper.alleleCounts[20]
	                              + helper.alleleCounts[65] + helper.alleleCounts[68];

	// If het-hets are 0 then do normal calculations
	// There is no phase uncertainty
	// Use phased math
#if TWK_SLAVE_DEBUG_MODE != 2
	if(number_of_hets == 0){
		helper.alleleCounts[TWK_LD_REFREF] = 2*helper.alleleCounts[0]  + helper.alleleCounts[1] + helper.alleleCounts[4] + helper.alleleCounts[16]  + helper.alleleCounts[64];
		helper.alleleCounts[TWK_LD_REFALT] = 2*helper.alleleCounts[5]  + helper.alleleCounts[4]  + helper.alleleCounts[1]  + helper.alleleCounts[21] + helper.alleleCounts[69];
		helper.alleleCounts[TWK_LD_ALTREF] = 2*helper.alleleCounts[80] + helper.alleleCounts[64] + helper.alleleCounts[16] + helper.alleleCounts[81] + helper.alleleCounts[84];
		helper.alleleCounts[TWK_LD_ALTALT] = 2*helper.alleleCounts[85] + helper.alleleCounts[84] + helper.alleleCounts[81] + helper.alleleCounts[69] + helper.alleleCounts[21];

		// Update counter
		//++this->no_uncertainty;

		// Reset
		cur_rcd.ChiSqModel = 0;

		// Use standard math
		return(this->PhasedMath(b1,p1,b2,p2));
	}
#endif

	/*double n1111 = helper.alleleCounts[0]; // 0
	double n1112 = helper.alleleCounts[1] + helper.alleleCounts[4]; // 1+4
	double n1122 = helper.alleleCounts[5]; // 5
	double n1211 = helper.alleleCounts[16] + helper.alleleCounts[64]; // 16 + 64
	double n1212 = helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68]; // 17 + 20 + 65 + 68
	double n1222 = helper.alleleCounts[21] + helper.alleleCounts[69]; // 21 + 69
	double n2211 = helper.alleleCounts[80]; // 80
	double n2212 = helper.alleleCounts[81] + helper.alleleCounts[84]; // 81 + 84
	double n2222 = helper.alleleCounts[85]; // 85
	*/
	//double n = (n1111 + n1112 + n1122 + n1211 + n1212 + n1222 + n2211 + n2212 + n2222);
	//double n = helper.totalHaplotypeCounts;
	const double P = ((helper.alleleCounts[0]+helper.alleleCounts[1] + helper.alleleCounts[4]+helper.alleleCounts[5])*2.0+(helper.alleleCounts[16] + helper.alleleCounts[64]+helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68]+helper.alleleCounts[21] + helper.alleleCounts[69]))/(2.0 * helper.totalHaplotypeCounts);
	const double Q = ((helper.alleleCounts[0]+helper.alleleCounts[16] + helper.alleleCounts[64]+helper.alleleCounts[80])*2.0+(helper.alleleCounts[1] + helper.alleleCounts[4]+helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68]+helper.alleleCounts[81] + helper.alleleCounts[84]))/(2.0 * helper.totalHaplotypeCounts);
	const double n11 = (2.0*helper.alleleCounts[0] + helper.alleleCounts[1] + helper.alleleCounts[4] + helper.alleleCounts[16] + helper.alleleCounts[64]);
	//const double n12 = (2.0*helper.alleleCounts[5] + helper.alleleCounts[1] + helper.alleleCounts[4] + helper.alleleCounts[21] + helper.alleleCounts[69]);
	//const double n21 = (2.0*helper.alleleCounts[80] + helper.alleleCounts[81] + helper.alleleCounts[84] + helper.alleleCounts[16] + helper.alleleCounts[64]);
	//const double n22 = (2.0*helper.alleleCounts[85] + helper.alleleCounts[81] + helper.alleleCounts[84] + helper.alleleCounts[21] + helper.alleleCounts[69]);
	const double minhap = n11 / (2.0 * helper.totalHaplotypeCounts);
	const double maxhap = (n11 + helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68]) / (2.0 * helper.totalHaplotypeCounts);
	//std::cerr << "haps=" << n11 << "," << n << "," << n1212 << std::endl;
	const double dee = -n11*P*Q;
	const double c = -n11*(1.0 - 2.0*P - 2.0*Q) - (helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68])*(1.0 - P - Q) + (2.0*helper.totalHaplotypeCounts*P*Q);
	const double b = 2.0*helper.totalHaplotypeCounts*(1.0 - 2.0*P - 2.0*Q) - 2.0*n11 - (helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68]);
	const double a = 4.0 * helper.totalHaplotypeCounts;

	/*std::cerr << n1111 << "\t" << n1112 << "\t" << n1122 << "\n"
			  << n1211 << "\t" << n1212 << "\t" << n1222 << "\n"
			  << n2211 << "\t" << n2212 << "\t" << n2222 << std::endl;
	std::cerr << "vals = " << n << "," << P << "," << Q << " haps= " << minhap << "-" << maxhap << "\t" << n11 << "," << n12 << "," << n21 << "," << n22 << std::endl;
	std::cerr << "vals2= " << a0 << "," << a1 << "," << a2 << "," << a3 << std::endl;
*/
	/*double a   = a3;
	double b   = a2;
	double c   = a1;
	double dee = a0;*/

	const double xN  = -b/(3.0*a);
	const double d2  = (pow(b,2)-3.0*a*c) / (9*pow(a,2));
	const double yN  = a * pow(xN,3) + b * pow(xN,2) + c * xN + dee;
	const double yN2 = pow(yN,2);
	const double h2  = 4 * pow(a,2) * pow(d2,3);

	/*const double p = ((helper.alleleCounts[TWK_LD_REFREF] + helper.alleleCounts[TWK_LD_REFALT]  + helper.alleleCounts[TWK_LD_ALTREF]  + helper.alleleCounts[TWK_LD_ALTALT])*2.0
				   + (helper.alleleCounts[16] + helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[21] + helper.alleleCounts[64] + helper.alleleCounts[65] + helper.alleleCounts[68] + helper.alleleCounts[69]))
				   / (2.0 * helper.totalHaplotypeCounts);
	const double q = ((helper.alleleCounts[TWK_LD_REFREF] + helper.alleleCounts[16] + helper.alleleCounts[64] + helper.alleleCounts[80])*2.0
				   + (helper.alleleCounts[TWK_LD_REFALT]  + helper.alleleCounts[TWK_LD_ALTREF]  + helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68] + helper.alleleCounts[81] + helper.alleleCounts[84]))
				   / (2.0 * helper.totalHaplotypeCounts);
	const double n11 = 2.0* helper.alleleCounts[TWK_LD_REFREF] + helper.alleleCounts[TWK_LD_REFALT] + helper.alleleCounts[TWK_LD_ALTREF] + helper.alleleCounts[16] + helper.alleleCounts[64];
*/
	// Not used for anything
	//const double n12 = (2.0*helper[5]  + helper[1]  + helper[4]  + helper[21] + helper[69]);
	//const double n21 = (2.0*helper[80] + helper[81] + helper[84] + helper[16] + helper[64]);
	//const double n22 = (2.0*helper[85] + helper[81] + helper[84] + helper[21] + helper[85]);

	/*////////////////////////
	// Cubic function: a3x^3 + a2x^2 + a1x + d = 0 <==> ax^3 + bx^2 + cx + d = 0
	// Cubic constants
	////////////////////////*/
	/*const double G   = 1.0 - 2.0*p - 2.0*q;
	const double dee = -n11*p*q;
	const double c   = -n11*G - number_of_hets*(1.0 - p - q) + 2.0*helper.totalHaplotypeCounts*p*q;
	const double b   = 2.0*helper.totalHaplotypeCounts*G - 2.0*n11 - number_of_hets;
	const double a   = 4.0 * helper.totalHaplotypeCounts;

	// Bounds for biological relevance
	const double minhap = n11 / (2.0 * helper.totalHaplotypeCounts);
	const double maxhap = (n11 + number_of_hets) / (2.0 * helper.totalHaplotypeCounts);

	// Cubic parameters
	const double xN  = -b / (3.0*a);
	const double d2  = (pow(b,2) - 3.0*a*c) / (9*pow(a,2));
	const double yN  = a * pow(xN,3) + b * pow(xN,2) + c * xN + dee;
	const double yN2 = pow(yN,2);
	const double h2  = 4 * pow(a,2) * pow(d2,3);*/

	// Difference between yN2 and h2
	const double __diff = yN2 - h2;

	/*std::cerr << "choose=" << yN2 << "," << h2 << " and " << abs(yN2-h2) << std::endl;
	if(abs(yN2-h2) <= 0.0000001){
		std::cerr << "no solution" << std::endl;
	}
	std::cerr << "pick=" << ((yN2 - h2) < 0) << ", " << ((yN2 - h2) > 0) << ", " << ((yN2 - h2) == 0) << std::endl;
*/
	// Begin cases
	if(__diff < 0) { // Yn2 < h2
		//option 3
		double h = pow(h2, 0.5);
		double theta = ((acos(-yN/h))/3.0);
		//delta = math.pow((yN/2.0*a),(1.0/3.0)) # is it correct to reuse this?
		double delta = pow(d2,0.5);
		double alpha = xN + 2.0 * delta * cos(theta);
		double beta  = xN + 2.0 * delta * cos(2.0 * M_PI/3.0 + theta);
		double gamma = xN + 2.0 * delta * cos(4.0 * M_PI/3.0 + theta);
		//result["p"] = p
		//result["q"] = q
		//result["pq"] = p * q

		//const double theta    = acos(-yN / sqrt(h2)) / 3.0;
		//const double constant = 2.0 * sqrt(d2);
		//const double alpha    = xN + constant * cos(theta);
		//const double beta     = xN + constant * cos(2.0*M_PI/3.0 + theta);
		//const double gamma    = xN + constant * cos(4.0*M_PI/3.0 + theta);

		uint8_t biologically_possible = 0;
		cur_rcd.ChiSqModel = std::numeric_limits<double>::max();
		double chosen = alpha;
		if(alpha >= minhap - TWK_ALLOWED_ROUNDING_ERROR && alpha <= maxhap + TWK_ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(alpha, P, Q);
			//chosen = &alpha;
			//std::cerr << "testing and chosing alpha: " << alpha << std::endl;
		}

		if(beta >= minhap - TWK_ALLOWED_ROUNDING_ERROR && beta <= maxhap + TWK_ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			//std::cerr << "testing beta" << std::endl;
			if(this->ChiSquaredUnphasedTable(beta, P, Q) < cur_rcd.ChiSqModel){
				chosen = beta;
				cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(beta, P, Q);
				//std::cerr << "chosen beta" << std::endl;
			}
		}

		if(gamma >= minhap - TWK_ALLOWED_ROUNDING_ERROR && gamma <= maxhap + TWK_ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			//std::cerr << "testing gamma" << std::endl;
			if(this->ChiSquaredUnphasedTable(gamma, P, Q) < cur_rcd.ChiSqModel){
				chosen = gamma;
				cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(gamma, P, Q); // implicit
				//std::cerr << "chosen gamma" << std::endl;
			}
		}

		if(biologically_possible == 0){
			//++this->impossible;
			return false;
		}

		if(biologically_possible > 1)
			cur_rcd.SetMultipleRoots();

		//++this->possible;
		return(this->ChooseF11Calculate(b1,p1,b2,p2,chosen, P, Q));

	} else if(__diff > 0){ // Yn2 > h2
		double number1 = 0.0;
		double number2 = 0.0;
		if((1.0/(2.0*a)*(-yN + pow((yN2 - h2),0.5))) < 0)
			 number1 = -pow(-(1.0/(2.0*a)*(-yN + pow((yN2 - h2),0.5))),1.0/3.0);
		else number1 = pow((1.0/(2.0*a)*(-yN + pow((yN2 - h2),0.5))),1.0/3.0);
		if((1.0/(2.0*a)*(-yN - pow((yN2 - h2),0.5))) < 0)
			 number2 = -pow(-(1.0/(2.0*a)*(-yN - pow((yN2 - h2),0.5))),1.0/3.0);
		else number2 = pow((1.0/(2.0*a)*(-yN - pow((yN2 - h2),0.5))),1.0/3.0);

		double alpha = xN + number1 + number2;

		if(!(alpha >= minhap - TWK_ALLOWED_ROUNDING_ERROR && alpha <= maxhap + TWK_ALLOWED_ROUNDING_ERROR)){
			//++this->impossible;
			return false;
		}

		//std::cerr << "chosen alpha: only available" << std::endl;
		cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(alpha, P, Q);
		//++this->possible;

		return(this->ChooseF11Calculate(b1,p1,b2,p2,alpha, P, Q));

	} else { // Yn2 == h2
		const double delta = pow((yN/2.0*a),(1.0/3.0));
		const double alpha = xN + delta; // alpha = beta in this case
		const double gamma = xN - 2.0*delta;

		if(std::isnan(alpha) || std::isnan(gamma)){
			//++this->impossible;
			return false;
		}

		uint8_t biologically_possible = 0;
		cur_rcd.ChiSqModel = std::numeric_limits<double>::max();
		double chosen = alpha;
		if(alpha >= minhap - TWK_ALLOWED_ROUNDING_ERROR && alpha <= maxhap + TWK_ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(alpha, P, Q);
			//chosen = &alpha;
			//std::cerr << "testing and chosing alpha/beta" << std::endl;
		}

		if(gamma >= minhap - TWK_ALLOWED_ROUNDING_ERROR && gamma <= maxhap + TWK_ALLOWED_ROUNDING_ERROR){
			++biologically_possible;
			//std::cerr << "testing gamma" << std::endl;
			if(this->ChiSquaredUnphasedTable(gamma, P, Q) < cur_rcd.ChiSqModel){
				chosen = gamma;
				cur_rcd.ChiSqModel = this->ChiSquaredUnphasedTable(gamma, P, Q); // implicit
				//std::cerr << "choisng gamma" << std::endl;
			}
		}

		if(biologically_possible == 0){
			//++this->impossible;
			return false;
		}

		//++this->possible;
		return(this->ChooseF11Calculate(b1,p1,b2,p2,chosen, P, Q));
	}
	return(false);
}

double twk_ld_engine::ChiSquaredUnphasedTable(const double target, const double p, const double q){
	const double f12   = p - target;
	const double f21   = q - target;
	const double f22   = 1 - (target + f12 + f21);

	const double e1111 = helper.totalHaplotypeCounts * pow(target,2);
	const double e1112 = 2 * helper.totalHaplotypeCounts * target * f12;
	const double e1122 = helper.totalHaplotypeCounts * pow(f12,2);
	const double e1211 = 2 * helper.totalHaplotypeCounts * target * f21;
	const double e1212 = 2 * helper.totalHaplotypeCounts * f12 * f21 + 2 * helper.totalHaplotypeCounts * target * f22;
	const double e1222 = 2 * helper.totalHaplotypeCounts * f12 * f22;
	const double e2211 = helper.totalHaplotypeCounts * pow(f21,2);
	const double e2212 = 2 * helper.totalHaplotypeCounts * f21 * f22;
	const double e2222 = helper.totalHaplotypeCounts * pow(f22,2);

	const double chisq1111 = e1111 > 0 ? pow((double)helper.alleleCounts[TWK_LD_REFREF]  - e1111, 2) / e1111 : 0,
				 chisq1112 = e1112 > 0 ? pow((double)helper.alleleCounts[TWK_LD_REFALT]  + helper.alleleCounts[TWK_LD_ALTREF] - e1112, 2) / e1112 : 0,
				 chisq1122 = e1122 > 0 ? pow((double)helper.alleleCounts[TWK_LD_ALTALT]  - e1122, 2) / e1122 : 0,
				 chisq1211 = e1211 > 0 ? pow((double)helper.alleleCounts[16] + helper.alleleCounts[64] - e1211, 2) / e1211 : 0,
				 chisq1212 = e1212 > 0 ? pow((double)helper.alleleCounts[17] + helper.alleleCounts[20] + helper.alleleCounts[65] + helper.alleleCounts[68] - e1212, 2) / e1212 : 0,
				 chisq1222 = e1222 > 0 ? pow((double)helper.alleleCounts[21] + helper.alleleCounts[69] - e1222, 2) / e1222 : 0,
				 chisq2211 = e2211 > 0 ? pow((double)helper.alleleCounts[80] - e2211, 2) / e2211 : 0,
				 chisq2212 = e2212 > 0 ? pow((double)helper.alleleCounts[81] + helper.alleleCounts[84] - e2212, 2) / e2212 : 0,
				 chisq2222 = e2222 > 0 ? pow((double)helper.alleleCounts[85] - e2222, 2) / e2222 : 0;

	return(chisq1111 + chisq1112 + chisq1122 + chisq1211 + chisq1212 + chisq1222 + chisq2211 + chisq2212 + chisq2222);
}

bool twk_ld_engine::ChooseF11Calculate(const twk1_ldd_blk& b1, const uint32_t pos1,
                                       const twk1_ldd_blk& b2, const uint32_t pos2,
                                       const double target,
                                       const double p, const double q)
{
	/*const double p1 = target;
	const double p2 = p - p1;
	const double q1 = q - p1;
	const double q2 = (1-(p1+p2+q1) < 0 ? 0 : 1-(p1+p2+q1));
*/
	double f11 = target;
	double f12 = p - f11;
	double f21 = q - f11;
	double f22 = 1 - (f11 + f12 + f21);
	double D = (f11 * f22) - (f12 * f21);
	double Dmax = 0;
	if(D >= 0.0) Dmax = std::min(p * (1.0-q), q * (1.0-p));
	else Dmax = std::min(p*q,(1-p)*(1-q));
	double Dprime = D / Dmax;
	double rsquared = (D * D) / (p * (1-p) * q * (1-q));

	cur_rcd.D  = D;
	cur_rcd.R2 = rsquared;

	//std::cerr << "unphased-choose=" << p1 << "," << p2 << "," << q1 << "," << q2 << " with p=" << p << " q=" << q << " R2=" << cur_rcd.R2  << std::endl;
	//if(rsquared > 0.5) exit(1);

	if(cur_rcd.R2 < settings.minR2 || cur_rcd.R2 > settings.maxR2){
		cur_rcd.controller = 0;
		return false;
	}

	cur_rcd.R  = sqrt(cur_rcd.R2);

	cur_rcd[TWK_LD_SIMD_REFREF] = f11 * 2*helper.totalHaplotypeCounts;
	cur_rcd[TWK_LD_SIMD_REFALT] = f12 * 2*helper.totalHaplotypeCounts;
	cur_rcd[TWK_LD_SIMD_ALTREF] = f21 * 2*helper.totalHaplotypeCounts;
	cur_rcd[TWK_LD_SIMD_ALTALT] = f22 * 2*helper.totalHaplotypeCounts;

	//std::cerr << p1 << "," << p2 << "," << q1 << "," << q2 << " and " << cur_rcd[TWK_LD_SIMD_REFREF] << "," << cur_rcd[TWK_LD_SIMD_REFALT] << "," << cur_rcd[TWK_LD_SIMD_ALTREF] << "," << cur_rcd[TWK_LD_SIMD_ALTALT] << " with R2=" << cur_rcd.R2 << std::endl;

	if(cur_rcd[TWK_LD_SIMD_REFREF] < cur_rcd[TWK_LD_SIMD_ALTALT]){
		if(cur_rcd[TWK_LD_SIMD_REFALT] + cur_rcd[TWK_LD_SIMD_ALTREF] + cur_rcd[TWK_LD_SIMD_REFREF] < 5){
			//std::cerr << "return=" << cur_rcd[TWK_LD_SIMD_REFALT] + cur_rcd[TWK_LD_SIMD_ALTREF] + cur_rcd[TWK_LD_SIMD_ALTALT] << std::endl;
			cur_rcd.controller = 0;
			return false;
		}
	} else {
		if(cur_rcd[TWK_LD_SIMD_ALTALT] + cur_rcd[TWK_LD_SIMD_REFALT] + cur_rcd[TWK_LD_SIMD_ALTREF] < 5){
			//std::cerr << "return=" << cur_rcd[TWK_LD_SIMD_REFALT] + cur_rcd[TWK_LD_SIMD_ALTREF] + cur_rcd[TWK_LD_SIMD_ALTALT] << std::endl;
			cur_rcd.controller = 0;
			return false;
		}
	}

	double dmax = 0;
	if(cur_rcd.D >= 0) dmax = p*(1.0-q) < q*(1.0-p) ? p*(1.0-q) : q*(1.0-p);
	else dmax = p*q < (1-p)*(1-q) ? -p*q : -(1-p)*(1-q);
	cur_rcd.Dprime = cur_rcd.D / dmax;

	if(cur_rcd.Dprime < settings.minDprime || cur_rcd.Dprime > settings.maxDprime){
		cur_rcd.controller = 0;
		return false;
	}

	double left,right,both;
	kt_fisher_exact(round(cur_rcd[TWK_LD_SIMD_REFREF]),round(cur_rcd[TWK_LD_SIMD_REFALT]),
					round(cur_rcd[TWK_LD_SIMD_ALTREF]),round(cur_rcd[TWK_LD_SIMD_ALTALT]),
					&left,&right,&both);
	cur_rcd.P = both;

	if(cur_rcd.P > settings.minP){
		cur_rcd.controller = 0;
		return false;
	}

	cur_rcd.Apos = b1.blk->rcds[pos1].pos;
	cur_rcd.Bpos = b2.blk->rcds[pos2].pos;
	cur_rcd.ridA = b1.blk->rcds[pos1].rid;
	cur_rcd.ridB = b2.blk->rcds[pos2].rid;
	cur_rcd.ChiSqModel = 0;
	cur_rcd.ChiSqFisher = (cur_rcd[TWK_LD_SIMD_REFREF] + cur_rcd[TWK_LD_SIMD_REFALT] + cur_rcd[TWK_LD_SIMD_ALTREF] + cur_rcd[TWK_LD_SIMD_ALTALT]) * cur_rcd.R2;

	//std::cerr << "unphased -> D=" << helper.D << ",Dprime=" << helper.Dprime << ",R2=" << helper.R2 << ",R=" << helper.R << std::endl;
	cur_rcd.SetLowACA(b1.blk->rcds[pos1].ac < TWK_LOW_AC_THRESHOLD);
	cur_rcd.SetLowACB(b2.blk->rcds[pos2].ac < TWK_LOW_AC_THRESHOLD);
	cur_rcd.SetCompleteLD(cur_rcd[TWK_LD_SIMD_REFREF] < 1 || cur_rcd[TWK_LD_SIMD_REFALT] < 1 || cur_rcd[TWK_LD_SIMD_ALTREF] < 1 || cur_rcd[TWK_LD_SIMD_ALTALT] < 1);
	cur_rcd.SetPerfectLD(cur_rcd.R2 > 0.99);
	cur_rcd.SetHasMissingValuesA(b1.blk->rcds[pos1].an);
	cur_rcd.SetHasMissingValuesB(b2.blk->rcds[pos2].an);
	int32_t diff = (int32_t)b1.blk->rcds[pos1].pos - b2.blk->rcds[pos2].pos;
	cur_rcd.SetLongRange(abs(diff) > TWK_LONG_RANGE_THRESHOLD && (cur_rcd.ridA == cur_rcd.ridB));
	cur_rcd.SetSameContig(cur_rcd.ridA == cur_rcd.ridB);
	cur_rcd.SetInvalidHWEA(b1.blk->rcds[pos1].hwe < TWK_INVALID_HWE_THRESHOLD);
	cur_rcd.SetInvalidHWEB(b2.blk->rcds[pos2].hwe < TWK_INVALID_HWE_THRESHOLD);

	//if(helper.R2 > 0.5) exit(1);

#if TWK_SLAVE_DEBUG_MODE == 3
	if(cur_rcd.R2 > 0.1){
		twk_debug_pos2   = b1.blk->rcds[pos1].pos;
		twk_debug_pos2_2 = b2.blk->rcds[pos2].pos;
		if(twk_debug_pos1 == twk_debug_pos2 && twk_debug_pos1_2 == twk_debug_pos2_2){
			std::cout << "P\t" << cur_rcd << '\n';
			std::cout << "U\t" << cur_rcd << '\n';
		}
		//std::cout << "U\t" << b1.blk->rcds[pos1].pos << "\t" << b2.blk->rcds[pos2].pos << "\t" << helper.D << "\t" << helper.Dprime << "\t" << helper.R << "\t" << helper.R2 << '\n';
	}
	#endif

	// If the number of rcds written is equal to the flush limit then
	// compress and write output.
	if(blk_f.n == n_lim || irecF.rid != cur_rcd.ridA || irecR.rid != cur_rcd.ridB){
		if(this->CompressBlock() == false) return false;
		irecF.rid    = cur_rcd.ridA;
		irecF.ridB   = cur_rcd.ridB;
		irecF.minpos = cur_rcd.Apos;
		irecF.maxpos = cur_rcd.Apos;
		irecR.rid    = cur_rcd.ridB;
		irecR.ridB   = cur_rcd.ridA;
		irecR.minpos = cur_rcd.Bpos;
		irecR.maxpos = cur_rcd.Bpos;
	}
	// If ridB is mixed then set index value to -1.
	if(irecF.ridB != cur_rcd.ridB) irecF.ridB = -1;
	if(irecR.ridB != cur_rcd.ridA) irecR.ridB = -1;

	// Update index
	irecF.maxpos = cur_rcd.Apos;
	irecR.maxpos = cur_rcd.Bpos;

	// Add forward.
	blk_f += cur_rcd;
	// Swap tuple (ridA:posA) with (ridB:posB).
	uint32_t temp = cur_rcd.Apos;
	cur_rcd.Apos = cur_rcd.Bpos;
	cur_rcd.Bpos = temp;
	std::swap(cur_rcd.ridA,cur_rcd.ridB);
	// Add reverse.
	blk_r += cur_rcd;
	// Update tickers.
	n_out += 2; //n_out_tick += 2;

	// Update ticker.
	//if(n_out_tick == 300){ progress->n_out += 300; n_out_tick = 0; }

	// Reset controller.
	cur_rcd.controller = 0;

	return false;
}

bool twk_ld_engine::CompressFwd(){
	if(blk_f.n){
		//progress->n_out += blk_f.n;
		ibuf << blk_f;

		if(zcodec.Compress(ibuf, obuf, settings.c_level) == false){
			std::cerr << "failed compression" << std::endl;
			return false;
		}

		t_out += ibuf.size() / sizeof(twk1_two_t);
		progress->b_out += ibuf.size();
		//std::cerr << "F=" << buf_f.size() << "->" << obuf.size() << " -> " << (float)buf_f.size()/obuf.size() << std::endl;
		//writer->flush(); // flush to make sure offset is good.
		//irecF.foff = writer->stream.tellp();
		irecF.b_cmp = obuf.size(); // keep here because writer resets obuf.
		writer->Add(ibuf.size(), obuf.size(), obuf, irecF);
		//writer->flush(); // flush to make sure offset is good.
		//irecF.fend = writer->stream.tellp();
		irecF.n = blk_f.n;
		irecF.b_unc = twk1_two_t::packed_size * irecF.n + 2*sizeof(uint32_t);
		index->AddThreadSafe(irecF);
		//std::cerr << irecF.n << "," << irecF.foff << "," << irecF.fend << "," << irecF.rid << ":" << irecF.minpos << "-" << irecF.maxpos << "," << irecF.ridB << std::endl;
		ibuf.reset();
		blk_f.reset();
		irecF.clear();
	}

	return true;
}

bool twk_ld_engine::CompressRev(){
	if(blk_r.n){
		ibuf << blk_r;
		//progress->n_out += blk_r.n;

		if(zcodec.Compress(ibuf, obuf, settings.c_level) == false){
			std::cerr << "failed compression" << std::endl;
			return false;
		}

		t_out += ibuf.size() / sizeof(twk1_two_t);
		progress->b_out += ibuf.size();
		//std::cerr << "F=" << buf_f.size() << "->" << obuf.size() << " -> " << (float)buf_f.size()/obuf.size() << std::endl;
		//writer->flush(); // flush to make sure offset is good.
		//irecR.foff = writer->stream.tellp();
		irecR.b_cmp = obuf.size(); // keep here because writer resets obuf.
		writer->Add(ibuf.size(), obuf.size(), obuf, irecR);
		//writer->flush(); // flush to make sure offset is good.
		//irecR.fend = writer->stream.tellp();
		irecR.n = blk_r.n;
		irecR.b_unc = twk1_two_t::packed_size * irecR.n + 2*sizeof(uint32_t);
		index->AddThreadSafe(irecR);
		//std::cerr << irecF.n << "," << irecF.foff << "," << irecF.fend << "," << irecF.rid << ":" << irecF.minpos << "-" << irecF.maxpos << "," << irecF.ridB << std::endl;
		ibuf.reset();
		blk_r.reset();
		irecR.clear();
	}

	return true;
}

bool twk_ld_engine::CompressBlock(){
	if(CompressFwd() == false) return false;
	if(CompressRev() == false) return false;
	//n_out = 0;

	return true;
}

/****************************
*  LD slave
****************************/
twk_ld_slave::twk_ld_slave() : n_s(0), n_total(0),
	i_start(0), j_start(0), prev_i(0), prev_j(0), n_cycles(0),
	ticker(nullptr), thread(nullptr), ldd(nullptr),
	progress(nullptr), settings(nullptr)
{}

twk_ld_slave::~twk_ld_slave(){ delete thread; }

std::thread* twk_ld_slave::Start(){
	delete thread; thread = nullptr;

	// Special case for single
	if(settings->single){
		thread = new std::thread(&twk_ld_slave::CalculateSingle, this, nullptr);
		return(thread);
	}

	if(settings->force_phased && settings->low_memory && settings->bitmaps)
		if(settings->window) thread = new std::thread(&twk_ld_slave::CalculatePhasedBitmapWindow, this, nullptr);
		else thread = new std::thread(&twk_ld_slave::CalculatePhasedBitmap, this, nullptr);
	else if(settings->force_phased){
		if(settings->window) thread = new std::thread(&twk_ld_slave::CalculatePhasedWindow, this, nullptr);
		else thread = new std::thread(&twk_ld_slave::CalculatePhased, this, nullptr);
	} else if(settings->forced_unphased){
		if(settings->window) thread = new std::thread(&twk_ld_slave::CalculateUnphasedWindow, this, nullptr);
		else thread = new std::thread(&twk_ld_slave::CalculateUnphased, this, nullptr);
	}
	else {
		thread = new std::thread(&twk_ld_slave::Calculate, this, nullptr);
	}

	return(thread);
}

void twk_ld_slave::UpdateBlocks(twk1_ldd_blk* blks,
                                const uint32_t& from,
                                const uint32_t& to)
{
	if(settings->low_memory == false)
		this->UpdateBlocksPreloaded(blks,from,to);
	else
		this->UpdateBlocksGenerate(blks,from,to);
}

void twk_ld_slave::UpdateBlocksPreloaded(twk1_ldd_blk* blks,
                                         const uint32_t& from,
                                         const uint32_t& to)
{
	const uint32_t add = ticker->diag ? 0 : (ticker->tL - ticker->fL);
	blks[0].SetPreloaded(ldd[from - i_start]);
	blks[1].SetPreloaded(ldd[add + (to - j_start)]);
	prev_i = from - i_start;
	prev_j = add + (to - j_start);
	++n_cycles;
}

void twk_ld_slave::UpdateBlocksGenerate(twk1_ldd_blk* blks,
                                        const uint32_t& from,
                                        const uint32_t& to)
{
	const uint32_t add = ticker->diag ? 0 : (ticker->tL - ticker->fL);

	if(n_cycles == 0 || prev_i != from - i_start){
		blks[0] = ldd[from - i_start];
		blks[0].Inflate(n_s, settings->ldd_load_type, true);
	} else {
		blks[0].blk   = ldd[from - i_start].blk;
		blks[0].n_rec = ldd[from - i_start].blk->n;
	}

	if(n_cycles == 0 || prev_j != add + (to - j_start)){
		blks[1] = ldd[add + (to - j_start)];
		blks[1].Inflate(n_s, settings->ldd_load_type, true);
	} else {
		blks[1].blk   = ldd[add + (to - j_start)].blk;
		blks[1].n_rec = ldd[add + (to - j_start)].blk->n;
	}

	prev_i = from - i_start;
	prev_j = add + (to - j_start);
	++n_cycles;
}

void twk_ld_slave::Phased(const twk1_t* rcds0,
                          const twk1_t* rcds1,
                          const twk1_ldd_blk* blocks,
                          const uint8_t type,
                          twk_ld_perf* perf)
{
	// Heuristically determined linear model at varying number of samples
	// using SSE4.2
	// y = 0.008145*n_s + 32.8227
	//const uint32_t thresh_nomiss = 0.008145*n_s + 32.8227;
	// RLEP-BVP intersection
	// y = 0.0047*n_s + 5.2913
	const uint32_t thresh_miss = 0.0047*n_s + 5.2913;
	//std::cerr << thresh_nomiss << "," << thresh_miss << std::endl;

	if(type == 1){
		uint64_t cur_out = engine.n_out;
		for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
			cur_out = engine.n_out;
			for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
				if(rcds0[i].ac + rcds0[j].ac <= 2){
					continue;
				}

				if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
					engine.PhasedListVector(blocks[0],i,blocks[0],j,perf);
				} else {
					if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
						engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
					else
						engine.PhasedVectorized(blocks[0],i,blocks[0],j,perf);
				}
			}
			progress->n_out += engine.n_out - cur_out;
		}
		progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
	} else {
		// Cache blocking
		uint32_t bsize = (256e3/2) / (2*n_s/8);
		bsize = (bsize == 0 ? 10 : bsize);
		const uint32_t n_blocks1 = blocks[0].n_rec / bsize;
		const uint32_t n_blocks2 = blocks[1].n_rec / bsize;

		uint64_t d = 0;
		uint64_t cur_out = engine.n_out;
		for(uint32_t ii = 0; ii < n_blocks1*bsize; ii += bsize){
			for(uint32_t jj = 0; jj < n_blocks2*bsize; jj += bsize){
				for(uint32_t i = ii; i < ii + bsize; ++i){
					cur_out = engine.n_out;
					for(uint32_t j = jj; j < jj + bsize; ++j){
						++d;
						if( rcds0[i].ac + rcds1[j].ac <= 2){
							continue;
						}

						//std::cerr << ii << "/" << n_blocks1 << "," << jj << "/" << n_blocks2 << "," << i << "/" << ii+bsize << "," << j << "/" << jj+bsize << " bsize=" << bsize << std::endl;
						if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
							engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
							else
								engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
						}
					}
					progress->n_out += engine.n_out - cur_out;
					//std::cerr << "end j" << std::endl;
				}
				//std::cerr << "end i" << std::endl;
			}
			// residual j that does not fit in a block
			for(uint32_t i = ii; i < ii + bsize; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = n_blocks2*bsize; j < blocks[1].n_rec; ++j){
					++d;
					if( rcds0[i].ac + rcds1[j].ac <= 2){
						continue;
					}

					if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
						engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
						else
							engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
					}

				}
				progress->n_out += engine.n_out - cur_out;
			}
		}

		//std::cerr << "i=" << n_blocks1*bsize << "/" << blocks[0].n_rec << std::endl;
		//std::cerr << "j=" << n_blocks2*bsize << "/" << blocks[1].n_rec << std::endl;

		for(uint32_t i = n_blocks1*bsize; i < blocks[0].n_rec; ++i){
			cur_out = engine.n_out;
			for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
				++d;
				if( rcds0[i].ac + rcds1[j].ac <= 2){
					continue;
				}

				if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
					engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
				} else {
					if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
						engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
					else
						engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
				}
			}
			progress->n_out += engine.n_out - cur_out;
		}

		//std::cerr << "done=" << d << "/" << blocks[0].n_rec * blocks[1].n_rec << std::endl;
		progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
	}
}

void twk_ld_slave::Unphased(const twk1_t* rcds0, const twk1_t* rcds1, const twk1_ldd_blk* blocks, const uint8_t type, twk_ld_perf* perf){
	// Heuristically determined linear model at varying number of samples
	// using SSE4.2
	// y = 0.0088*n_s + 18.972
	const uint32_t thresh_nomiss = 0.0088*n_s + 18.972;
	// RLEU-BVU intersection
	// y = 0.012*n_s + 22.3661
	const uint32_t thresh_miss = 0.012*n_s + 22.3661;

	if(type == 1){
		uint64_t cur_out = engine.n_out;
		for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
			cur_out = engine.n_out;
			for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
				if(rcds0[i].ac + rcds0[j].ac <= 2){
					continue;
				}

#if(TWK_SLAVE_DEBUG_MODE == 2)
				if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
					engine.UnphasedRunlength(blocks[0],i,blocks[0],j,perf);
					engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[0],j,perf);
				} else {
					engine.UnphasedRunlength(blocks[0],i,blocks[0],j,perf);
					engine.UnphasedVectorized(blocks[0],i,blocks[0],j,perf);
				}
#else
				if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
					if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
						engine.UnphasedRunlength(blocks[0],i,blocks[0],j,perf);
					else {
						engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[0],j,perf);
					}
				} else {
					if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
						engine.UnphasedRunlength(blocks[0],i,blocks[0],j,perf);
					else
						engine.UnphasedVectorized(blocks[0],i,blocks[0],j,perf);
				}
#endif
			}
			progress->n_out += engine.n_out - cur_out;
		}
		progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
	} else {
		uint32_t bsize = (256e3/2) / (2*n_s/8);
		bsize = (bsize == 0 ? 10 : bsize);
		const uint32_t n_blocks1 = blocks[0].n_rec / bsize;
		const uint32_t n_blocks2 = blocks[1].n_rec / bsize;

		uint64_t d = 0;
		uint64_t cur_out = engine.n_out;
		for(uint32_t ii = 0; ii < n_blocks1*bsize; ii += bsize){
			for(uint32_t jj = 0; jj < n_blocks2*bsize; jj += bsize){
				for(uint32_t i = ii; i < ii + bsize; ++i){
					cur_out = engine.n_out;
					for(uint32_t j = jj; j < jj + bsize; ++j){
						++d;
						if( rcds0[i].ac + rcds1[j].ac <= 2){
							continue;
						}

						//std::cerr << ii << "/" << n_blocks1 << "," << jj << "/" << n_blocks2 << "," << i << "/" << ii+bsize << "," << j << "/" << jj+bsize << " bsize=" << bsize << std::endl;
#if(TWK_SLAVE_DEBUG_MODE == 2)
						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
							engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
						} else {
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
							engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
						}
#else
						if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
							if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
								engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
							else {
								engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
							}
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
							else
								engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
						}
#endif
					}
					progress->n_out += engine.n_out - cur_out;
					//std::cerr << "end j" << std::endl;
				}
				//std::cerr << "end i" << std::endl;
			}
			// residual j that does not fit in a block
			for(uint32_t i = ii; i < ii + bsize; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = n_blocks2*bsize; j < blocks[1].n_rec; ++j){
					++d;
					if( rcds0[i].ac + rcds1[j].ac <= 2){
						continue;
					}

#if(TWK_SLAVE_DEBUG_MODE == 2)
					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
						engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
					} else {
						engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
						engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
					}
#else
					if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
						if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
						else {
							engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
						}
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
						else
							engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
					}
#endif

				}
				progress->n_out += engine.n_out - cur_out;
			}
		}

		//std::cerr << "i=" << n_blocks1*bsize << "/" << blocks[0].n_rec << std::endl;
		//std::cerr << "j=" << n_blocks2*bsize << "/" << blocks[1].n_rec << std::endl;

		for(uint32_t i = n_blocks1*bsize; i < blocks[0].n_rec; ++i){
			cur_out = engine.n_out;
			for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
				++d;
				if( rcds0[i].ac + rcds1[j].ac <= 2){
					continue;
				}

#if(TWK_SLAVE_DEBUG_MODE == 2)
				if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
					engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
					engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
				} else {
					engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
					engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
				}
#else
				if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
					if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
						engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
					else {
						engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
					}
				} else {
					if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
						engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
					else
						engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
				}
#endif

			}
			progress->n_out += engine.n_out - cur_out;
		}

		//std::cerr << "done=" << d << "/" << blocks[0].n_rec * blocks[1].n_rec << std::endl;
		progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
	}
}

bool twk_ld_slave::CalculatePhased(twk_ld_perf* perf){
	twk1_ldd_blk blocks[2];
	uint32_t from, to; uint8_t type;
	Timer timer; timer.Start();

	i_start = ticker->fL; j_start = ticker->fR;
	prev_i = 0; prev_j = 0;
	n_cycles = 0;
	const twk1_t* rcds0 = nullptr;
	const twk1_t* rcds1 = nullptr;

	while(true){
		if(!ticker->Get(from, to, type)) break;
		this->UpdateBlocks(blocks,from,to);

		rcds0 = blocks[0].blk->rcds;
		rcds1 = blocks[1].blk->rcds;

		// Compute phased math
		Phased(rcds0, rcds1, blocks, type, perf);
	}
	//std::cerr << "done" << std::endl;

	// if preloaded
	if(settings->low_memory == false){
		blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
		blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
	}

	// compress last
	if(this->engine.CompressBlock() == false)
		return false;

	return true;
}

bool twk_ld_slave::CalculateSingle(twk_ld_perf* perf){
	twk1_ldd_blk blocks[2];
	uint32_t from, to; uint8_t type;
	Timer timer; timer.Start();

	i_start = ticker->fL; j_start = ticker->fR;
	prev_i = 0; prev_j = 0;
	n_cycles = 0;

	// Heuristically determined
	uint32_t cycle_thresh_p = n_s / 45 == 0 ? 2 : n_s / 45;
	if(settings->cycle_threshold != 0)
		cycle_thresh_p = settings->cycle_threshold;
	uint32_t cycle_thresh_u = (n_s / 60 == 0 ? 2 : n_s / 60);
	if(settings->cycle_threshold != 0)
		cycle_thresh_u = settings->cycle_threshold;

	const uint32_t thresh_miss = 0.0047*n_s + 5.2913;

	const twk1_t* rcds0 = nullptr;
	const twk1_t* rcds1 = nullptr;

	while(true){
		if(!ticker->Get(from, to, type)) break;
		//this->UpdateBlocks(blocks,from,to);
		blocks[0].SetPreloaded(ldd[0]);
		blocks[1].SetPreloaded(ldd[to]);
		if(to == 0) type = 1;
		else type = 0;

		//std::cerr << "f/t/t: " << from << "," << to << "," << (int)type << " of " << ticker->fL << "-" << ticker->tR << std::endl;
		//std::cerr << "starts=" << i_start << "/" << j_start << std::endl;
		//std::cerr << "sizes=" << blocks[0].blk->n << "," << blocks[1].blk->n << std::endl;
		//std::cerr << "range=" << blocks[0].blk->minpos << "-" << blocks[0].blk->maxpos << " and " << blocks[1].blk->minpos << "-" << blocks[1].blk->maxpos << std::endl;

		rcds0 = blocks[0].blk->rcds;
		rcds1 = blocks[1].blk->rcds;

		if(type == 1){
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
					/*if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac <= 2){
						continue;
					}*/

					if(blocks[0].blk->rcds[i].an || blocks[0].blk->rcds[j].an){
						if(blocks[0].blk->rcds[i].gt->n + blocks[0].blk->rcds[j].gt->n < cycle_thresh_u){
							engine.UnphasedRunlength(blocks[0],i,blocks[0],j,nullptr);
						} else {
							engine.UnphasedVectorized(blocks[0],i,blocks[0],j,nullptr);
						}
					} else {
						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedListVector(blocks[0],i,blocks[0],j,perf);
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
							else
								engine.PhasedVectorized(blocks[0],i,blocks[0],j,perf);
						}
					}
				}
			}
			progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
		} else {
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
					/*if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac <= 2 ){
						continue;
					}*/

					if(blocks[0].blk->rcds[i].an || blocks[1].blk->rcds[j].an){
						if(blocks[0].blk->rcds[i].gt->n + blocks[1].blk->rcds[j].gt->n < cycle_thresh_u)
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,nullptr);
						else {
							engine.UnphasedVectorized(blocks[0],i,blocks[1],j,nullptr);
						}
					} else {
						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
							else
								engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
						}
					}

				}
			}
			progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
		}
	}
	//std::cerr << "done" << std::endl;

	// if preloaded
	if(settings->low_memory == false){
		blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
		blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
	}

	// compress last
	if(this->engine.CompressBlock() == false)
		return false;

	return true;
}

bool twk_ld_slave::CalculateUnphased(twk_ld_perf* perf){
	twk1_ldd_blk blocks[2];
	uint32_t from, to; uint8_t type;
	Timer timer; timer.Start();

	i_start = ticker->fL; j_start = ticker->fR;
	prev_i = 0; prev_j = 0;
	n_cycles = 0;
	const twk1_t* rcds0 = nullptr;
	const twk1_t* rcds1 = nullptr;

	while(true){
		if(!ticker->Get(from, to, type)) break;

		this->UpdateBlocks(blocks,from,to);
		rcds0 = blocks[0].blk->rcds;
		rcds1 = blocks[1].blk->rcds;

		// Compute unphased
		Unphased(rcds0, rcds1, blocks, type, perf);
	}
	//std::cerr << "done" << std::endl;

	// if preloaded
	if(settings->low_memory == false){
		blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
		blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
	}

	// compress last
	if(this->engine.CompressBlock() == false)
		return false;

	return true;
}

bool twk_ld_slave::CalculatePhasedBitmap(twk_ld_perf* perf){
	twk1_ldd_blk blocks[2];
	uint32_t from, to; uint8_t type;
	Timer timer; timer.Start();

	i_start = ticker->fL; j_start = ticker->fR;
	prev_i = 0; prev_j = 0;
	n_cycles = 0;

	const twk1_t* rcds0 = nullptr;
	const twk1_t* rcds1 = nullptr;


	while(true){
		if(!ticker->Get(from, to, type)) break;

		this->UpdateBlocks(blocks,from,to);
		rcds0 = blocks[0].blk->rcds;
		rcds1 = blocks[1].blk->rcds;

		uint32_t cur_out = engine.n_out;
		if(type == 1){
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
					if(rcds0[i].ac + rcds0[j].ac <= 2){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						engine.PhasedBitmap(blocks[0],i,blocks[0],j,perf);
					} else {
						engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
		} else {
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
					if( rcds0[i].ac + rcds1[j].ac <= 2){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						engine.PhasedBitmap(blocks[0],i,blocks[1],j,perf);
					} else {
						engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
		}
	}

	// if preloaded
	if(settings->low_memory == false){
		blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
		blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
	}

	// compress last
	if(this->engine.CompressBlock() == false)
		return false;

	return true;
}

bool twk_ld_slave::CalculatePhasedBitmapWindow(twk_ld_perf* perf){
	twk1_ldd_blk blocks[2];
	uint32_t from, to; uint8_t type;
	Timer timer; timer.Start();

	i_start = ticker->fL; j_start = ticker->fR;
	prev_i = 0; prev_j = 0;
	n_cycles = 0;

	const twk1_t* rcds0 = nullptr;
	const twk1_t* rcds1 = nullptr;


	while(true){
		if(!ticker->Get(from, to, type)) break;

		this->UpdateBlocks(blocks,from,to);
		rcds0 = blocks[0].blk->rcds;
		rcds1 = blocks[1].blk->rcds;

		uint32_t cur_out = engine.n_out;
		if(type == 1){
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
					if(blocks[0].blk->rcds[i].rid != blocks[0].blk->rcds[j].rid
					   && (blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
					{
						//std::cerr << "oor=" << blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[j].pos << std::endl;
						continue;
					}

					if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac <= 2){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						engine.PhasedBitmap(blocks[0],i,blocks[0],j,perf);
					} else {
						engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
		} else {
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
					if(blocks[0].blk->rcds[i].rid != blocks[1].blk->rcds[j].rid
					   && (blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
					{
						//std::cerr << "oor=" << blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
						continue;
					}

					if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac <= 2 ){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						engine.PhasedBitmap(blocks[0],i,blocks[1],j,perf);
					} else {
						engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
		}
	}

	// if preloaded
	if(settings->low_memory == false){
		blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
		blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
	}

	// compress last
	if(this->engine.CompressBlock() == false)
		return false;

	return true;
}

bool twk_ld_slave::CalculatePhasedWindow(twk_ld_perf* perf){
	twk1_ldd_blk blocks[2];
	uint32_t from, to; uint8_t type;
	Timer timer; timer.Start();

	i_start = ticker->fL; j_start = ticker->fR;
	prev_i = 0; prev_j = 0;
	n_cycles = 0;

	const uint32_t thresh_miss = 0.0047*n_s + 5.2913;

	const twk1_t* rcds0 = nullptr;
	const twk1_t* rcds1 = nullptr;

	while(true){
		if(!ticker->Get(from, to, type)) break;

		this->UpdateBlocks(blocks,from,to);
		rcds0 = blocks[0].blk->rcds;
		rcds1 = blocks[1].blk->rcds;

		uint32_t cur_out = engine.n_out;
		if(type == 1){
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
					//std::cerr << "here=" << i << "/" << j << " diff=" << blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
					if(blocks[0].blk->rcds[i].rid == blocks[0].blk->rcds[j].rid
					   && (blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
					{
						//std::cerr << "oor=" << blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
						//progress->n_var += i*j;
						goto end_cycle;
						//continue;
					}

					if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac <= 2){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						engine.PhasedListVector(blocks[0],i,blocks[0],j,perf);
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
						else
							engine.PhasedVectorized(blocks[0],i,blocks[0],j,perf);
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
		} else {
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
					if(blocks[0].blk->rcds[i].rid == blocks[1].blk->rcds[j].rid
					   && (blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
					{
						//std::cerr << "oor=" << blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
						goto end_cycle;
						//continue;
					}

					if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac <= 2 ){
						continue;
					}

					if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
						engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
						else
							engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
		}
		end_cycle:
		progress->n_var += 0;
	}

	// if preloaded
	if(settings->low_memory == false){
		blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
		blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
	}

	// compress last
	if(this->engine.CompressBlock() == false)
		return false;

	return true;
}

bool twk_ld_slave::CalculateUnphasedWindow(twk_ld_perf* perf){
	std::cerr << "in unphased window" << std::endl;
	twk1_ldd_blk blocks[2];
	uint32_t from, to; uint8_t type;
	Timer timer; timer.Start();

	i_start = ticker->fL; j_start = ticker->fR;
	prev_i = 0; prev_j = 0;
	n_cycles = 0;

	// Heuristically determined linear model at varying number of samples
	// using SSE4.2
	// y = 0.0088*n_s + 18.972
	const uint32_t thresh_nomiss = 0.0088*n_s + 18.972;
	// RLEU-BVU intersection
	// y = 0.012*n_s + 22.3661
	const uint32_t thresh_miss = 0.012*n_s + 22.3661;

	const twk1_t* rcds0 = nullptr;
	const twk1_t* rcds1 = nullptr;

	while(true){
		if(!ticker->Get(from, to, type)) break;

		this->UpdateBlocks(blocks,from,to);
		rcds0 = blocks[0].blk->rcds;
		rcds1 = blocks[1].blk->rcds;

		uint32_t cur_out = engine.n_out;
		if(type == 1){
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
					//std::cerr << "here=" << i << "/" << j << " diff=" << blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
					if(blocks[0].blk->rcds[i].rid == blocks[0].blk->rcds[j].rid
					   && (blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
					{
						//std::cerr << "oor=" << blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
						//progress->n_var += i*j;
						goto end_cycle;
						//continue;
					}

					if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac <= 2){
						continue;
					}

					if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
						if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
							engine.UnphasedRunlength(blocks[0],i,blocks[0],j,perf);
						else {
							engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[0],j,perf);
						}
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.UnphasedRunlength(blocks[0],i,blocks[0],j,perf);
						else
							engine.UnphasedVectorized(blocks[0],i,blocks[0],j,perf);
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
		} else {
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
					if(blocks[0].blk->rcds[i].rid == blocks[1].blk->rcds[j].rid
					   && (blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) > settings->l_window)
					{
						//std::cerr << "oor=" << blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos << std::endl;
						goto end_cycle;
						//continue;
					}

					if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac <= 2 ){
						continue;
					}

					if(rcds0[i].gt_missing == false && rcds1[j].gt_missing == false){
						if(std::min(rcds0[i].ac, rcds0[j].ac) < thresh_nomiss)
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
						else {
							engine.UnphasedVectorizedNoMissing(blocks[0],i,blocks[1],j,perf);
						}
					} else {
						if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,perf);
						else
							engine.UnphasedVectorized(blocks[0],i,blocks[1],j,perf);
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
		}
		end_cycle:
		progress->n_var += 0;
	}

	// if preloaded
	if(settings->low_memory == false){
		blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
		blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
	}

	// compress last
	if(this->engine.CompressBlock() == false)
		return false;

	return true;
}

bool twk_ld_slave::Calculate(twk_ld_perf* perf){
	twk1_ldd_blk blocks[2];
	uint32_t from, to; uint8_t type;
	Timer timer; timer.Start();

	i_start = ticker->fL; j_start = ticker->fR;
	prev_i = 0; prev_j = 0;
	n_cycles = 0;

	// Heuristically determined
	uint32_t cycle_thresh_p = n_s / 45 == 0 ? 2 : n_s / 45;
	if(settings->cycle_threshold != 0)
		cycle_thresh_p = settings->cycle_threshold;
	uint32_t cycle_thresh_u = (n_s / 60 == 0 ? 2 : n_s / 60);
	if(settings->cycle_threshold != 0)
		cycle_thresh_u = settings->cycle_threshold;

	const uint32_t thresh_miss = 0.0047*n_s + 5.2913;

	const twk1_t* rcds0 = nullptr;
	const twk1_t* rcds1 = nullptr;

	while(true){
		if(!ticker->Get(from, to, type)) break;

		this->UpdateBlocks(blocks,from,to);
		rcds0 = blocks[0].blk->rcds;
		rcds1 = blocks[1].blk->rcds;

		uint32_t cur_out = engine.n_out;
		if(type == 1){
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
					if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac <= 2){
						continue;
					}

					if(blocks[0].blk->rcds[i].an || blocks[0].blk->rcds[j].an){
						if(blocks[0].blk->rcds[i].gt->n + blocks[0].blk->rcds[j].gt->n < cycle_thresh_u){
							engine.UnphasedRunlength(blocks[0],i,blocks[0],j,nullptr);
						} else {
							engine.UnphasedVectorized(blocks[0],i,blocks[0],j,nullptr);
						}
					} else {
						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedListVector(blocks[0],i,blocks[0],j,perf);
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.PhasedRunlength(blocks[0],i,blocks[0],j,perf);
							else
								engine.PhasedVectorized(blocks[0],i,blocks[0],j,perf);
						}
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
		} else {
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				cur_out = engine.n_out;
				for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
					if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac <= 2 ){
						continue;
					}

					if(blocks[0].blk->rcds[i].an || blocks[1].blk->rcds[j].an){
						if(blocks[0].blk->rcds[i].gt->n + blocks[1].blk->rcds[j].gt->n < cycle_thresh_u)
							engine.UnphasedRunlength(blocks[0],i,blocks[1],j,nullptr);
						else {
							engine.UnphasedVectorized(blocks[0],i,blocks[1],j,nullptr);
						}
					} else {
						if((rcds0[i].gt_missing || rcds1[j].gt_missing) == false){
							engine.PhasedListVector(blocks[0],i,blocks[1],j,perf);
						} else {
							if(rcds0[i].ac + rcds1[j].ac < thresh_miss)
								engine.PhasedRunlength(blocks[0],i,blocks[1],j,perf);
							else
								engine.PhasedVectorized(blocks[0],i,blocks[1],j,perf);
						}
					}
				}
				progress->n_out += engine.n_out - cur_out;
			}
			progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
		}

	}

	// if preloaded
	if(settings->low_memory == false){
		blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
		blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
	}

	// compress last
	if(this->engine.CompressBlock() == false)
		return false;

	return true;
}

bool twk_ld_slave::CalculatePerformance(twk_ld_engine::func f, twk_ld_perf* perf){
	twk1_ldd_blk blocks[2];
	uint32_t from, to; uint8_t type;
	Timer timer; timer.Start();

	i_start = ticker->fL; j_start = ticker->fR;
	prev_i = 0; prev_j = 0;
	n_cycles = 0;
	const twk1_t* rcds0 = nullptr;
	const twk1_t* rcds1 = nullptr;

	while(true){
		if(!ticker->Get(from, to, type)) break;

		this->UpdateBlocks(blocks,from,to);
		rcds0 = blocks[0].blk->rcds;
		rcds1 = blocks[1].blk->rcds;

		if(type == 1){
			for(uint32_t i = 0; i < blocks[0].n_rec; ++i){
				for(uint32_t j = i+1; j < blocks[0].n_rec; ++j){
					(engine.*f)(blocks[0],i,blocks[0],j,perf);
				}
			}
			progress->n_var += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
		} else {
			uint32_t bsize = (256e3/2) / (2*n_s/8);
			bsize = bsize < 1 ? 1 : bsize; // minimum block size
			const uint32_t n_blocks1 = blocks[0].n_rec / bsize;
			const uint32_t n_blocks2 = blocks[1].n_rec / bsize;

			// Cache blocking
			for(uint32_t ii = 0; ii < n_blocks1*bsize; ii += bsize){
				for(uint32_t jj = 0; jj < n_blocks2*bsize; jj += bsize){
					for(uint32_t i = ii; i < ii + bsize; ++i){
						for(uint32_t j = jj; j < jj + bsize; ++j){
							(engine.*f)(blocks[0],i,blocks[1],j,perf);
						}
					}
				}
				// residual j that does not fit in a block
				for(uint32_t i = ii; i < ii + bsize; ++i){
					for(uint32_t j = n_blocks2*bsize; j < blocks[1].n_rec; ++j){
						(engine.*f)(blocks[0],i,blocks[1],j,perf);
					}
				}
			}

			// Residual
			for(uint32_t i = n_blocks1*bsize; i < blocks[0].n_rec; ++i){
				for(uint32_t j = 0; j < blocks[1].n_rec; ++j){
					(engine.*f)(blocks[0],i,blocks[1],j,perf);
				}
			}
			progress->n_var += blocks[0].n_rec * blocks[1].n_rec;
		}
	}

	// if preloaded
	if(settings->low_memory == false){
		blocks[0].vec = nullptr; blocks[0].list = nullptr; blocks[0].bitmap = nullptr;
		blocks[1].vec = nullptr; blocks[1].list = nullptr; blocks[1].bitmap = nullptr;
	}

	return true;
}

}
