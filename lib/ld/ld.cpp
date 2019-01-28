#include <chrono>
#include <cstdlib>

#include "ld.h"
#include "intervals.h"
#include "fisher_math.h"
#include "ld/ld_engine.h"
#include "timer.h"

namespace tomahawk {

// Implementation to hide ldd from users. This approach allows
// us to not having to include third party headers in our releases.
class twk_ld::twk_ld_impl {
public:
	twk_ld_impl() : n_blks(0), m_blks(0), n_vnts(0), n_tree(0), ldd(nullptr), ldd2(nullptr){}
	~twk_ld_impl(){ delete[] ldd; delete[] ldd2; }

	/**<
	 * Wrapper function for parsing interval strings.
	 * @param reader Reference twk_reader instance.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseIntervalStrings(twk_reader& reader, const twk_ld_settings& settings){ return(intervals.ParseIntervalStrings(settings.ival_strings, reader.hdr)); }

	/**<
	 * Reads the desired tomahawk blocks given the balancer intervals. Internally
	 * decides if the block slicing is based on the universal set of blocks or a
	 * targetted subset (as decided by the interval slicing operation).
	 * @param reader
	 * @param bit
	 * @param balancer
	 * @param load
	 * @return
	 */
	bool LoadBlocks(twk_reader& reader,
					twk1_blk_iterator& bit,
					const twk_ld_balancer& balancer,
					const twk_ld_settings& settings);

	/**<
	 * Loading twk blocks for a single variant and its surrounding variants within
	 * some distance and on the same chromosome. This function loads the target variant
	 * in a single block and the other variants in separate blocks as usual. The
	 * identity of the target site is parameterized in the settings object.
	 * @param reader
	 * @param bit
	 * @param balancer
	 * @param load
	 * @return
	 */
	bool LoadTargetSingle(twk_reader& reader,
						  twk1_blk_iterator& bit,
						  const twk_ld_settings& settings,
						  const uint8_t load);

	/**<
	 * Construct interval container and trees given the pre-provided interval
	 * strings.
	 * @param reader
	 * @param bit
	 * @param load
	 * @return
	 */
	bool BuildIntervals(twk_reader& reader, const twk_ld_settings& settings, twk1_blk_iterator& bit);

	/**<
	 * Loads only the target blocks that overlap with the given vector of interval
	 * tuples as parameterized in the settings object.
	 * @param reader   Reference to twk reader.
	 * @param bit      Reference to a twk block iterator.
	 * @param balancer Reference of a pre-computed load balancer.
	 * @param settings Reference of a user-paramterized settings object.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool LoadTargetBlocks(twk_reader& reader,
						  twk1_blk_iterator& bit,
						  const twk_ld_balancer& balancer,
						  const twk_ld_settings& settings);

	/**<
	 * Loads all available twk blocks into memory. Internally spawns the maximum
	 * possible number of unpacking threads possible (as parameterized by settings).
	 * @param reader   Reference to twk reader.
	 * @param bit      Reference to a twk block iterator.
	 * @param balancer Reference of a pre-computed load balancer.
	 * @param settings Reference of a user-paramterized settings object.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool LoadAllBlocks(twk_reader& reader,
					   twk1_blk_iterator& bit,
					   const twk_ld_balancer& balancer,
					   const twk_ld_settings& settings);

public:
	uint32_t n_blks, m_blks, n_vnts, n_tree;
	twk_intervals intervals;
	twk1_ldd_blk* ldd;
	twk1_block_t* ldd2;
};


/****************************
*  twk_ld
****************************/
twk_ld::twk_ld() : mImpl(new twk_ld_impl)
{
}

twk_ld::~twk_ld(){
	delete mImpl;
}

bool twk_ld::twk_ld_impl::LoadBlocks(twk_reader& reader,
		twk1_blk_iterator& bit,
		const twk_ld_balancer& balancer,
		const twk_ld_settings& settings)
{
	if(intervals.overlap_blocks.size()) return(this->LoadTargetBlocks(reader, bit, balancer, settings));
	return(this->LoadAllBlocks(reader, bit, balancer, settings));
}

bool twk_ld::twk_ld_impl::LoadTargetSingle(twk_reader& reader,
		twk1_blk_iterator& bit,
		const twk_ld_settings& settings,
		const uint8_t load)
{
	if(settings.ival_strings.size() == 0){ // if have interval strings
		return false;
	}

	this->intervals.ivecs.resize(reader.hdr.GetNumberContigs());
	if(this->ParseIntervalStrings(reader, settings) == false)
		return false;

	uint32_t n_vecs = 0, ivec_rid = 0;
	for(int i = 0; i < intervals.ivecs.size(); ++i){
		//std::cerr << i << "/" << intervals.ivecs.size() << "\t" << intervals.ivecs[i].size() << std::endl;
		n_vecs += intervals.ivecs[i].size();
		if(intervals.ivecs[i].size()) ivec_rid = i;
	}
	assert(n_vecs == 1);
	//std::cerr << "target=" << intervals.ivecs[ivec_rid].front().start << "-" << intervals.ivecs[ivec_rid].front().stop << " val=" << intervals.ivecs[ivec_rid].front().value << std::endl;
	// Add flanking intervals
	int32_t f   = (int32_t)intervals.ivecs[ivec_rid][0].start - settings.l_surrounding;
	int32_t f_s = (int32_t)intervals.ivecs[ivec_rid][0].start - 1;
	f   = std::max(f,   0);
	f_s = std::max(f_s, 0);
	uint32_t g   = f; // stupid constructor design in intervalTree.h
	uint32_t g_s = f_s;
	uint32_t t = (int32_t)intervals.ivecs[ivec_rid][0].stop + settings.l_surrounding;
	intervals.ivecs[ivec_rid].push_back(twk_intervals::interval(g, g_s, 1)); // right side is inclusive here
	intervals.ivecs[ivec_rid].push_back(twk_intervals::interval(intervals.ivecs[ivec_rid][0].stop, t, 2));

	// Debug print.
	/*for(int i = 1; i < intervals.ivecs[ivec_rid].size(); ++i){
		std::cerr << "others=" << intervals.ivecs[ivec_rid][i].start << "-" << intervals.ivecs[ivec_rid][i].stop << " val=" << intervals.ivecs[ivec_rid][i].value << std::endl;
	}*/

	// Build intervals container.
	if(this->intervals.Build(reader.hdr.GetNumberContigs(), reader.index) == false){
		return false;
	}

	//ldd    = new twk1_ldd_blk[3];
	uint32_t ldd2_n = 1, ldd2_m = 1000;
	ldd2   = new twk1_block_t[ldd2_m];

	// First block has only 1 variant: the reference.
	// All other blocks
	/*std::cerr << "overlap blocks=" << intervals.overlap_blocks.size() << " and " << intervals.ivecs[ivec_rid].size() << std::endl;
	// Debug print.
	for(int i = 0; i < intervals.ivecs[ivec_rid].size(); ++i){
		std::cerr << "others=" << intervals.ivecs[ivec_rid][i].start << "-" << intervals.ivecs[ivec_rid][i].stop << " val=" << intervals.ivecs[ivec_rid][i].value << std::endl;
	}*/

	for(int i = 0; i < intervals.overlap_blocks.size(); ++i){
		// Make sure we don't seek to the same block twice. This will result
		// in incorrect duplicatation of variants.
		if(i!=0){
			if(intervals.overlap_blocks[i]->foff == intervals.overlap_blocks[i-1]->foff)
				continue;
		}

		bit.stream->seekg(intervals.overlap_blocks[i]->foff);
		if(bit.NextBlock() == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
			return false;
		}

		for(int j = 0; j < bit.blk.n; ++j){
			std::vector<twk_intervals::interval> mivals = intervals.itree[bit.blk.rcds[j].rid]->findOverlapping(bit.blk.rcds[j].pos+1, bit.blk.rcds[j].pos+1); // 1-base matching.
			// It is possible we get >1 overlapping interval matches since the
			// intervalTree checks for inclusive ranges [A,B]. In this case we
			// record the overlap as beloning to the 0-th bin or else is corrupted.
			if(mivals.size()){
				if(mivals.size() == 1){
					if(mivals[0].value == 0){
						if(ldd2[0].n) assert(ldd2[0].rcds[0].pos != bit.blk.rcds[j].pos);
						ldd2[0].Add(bit.blk.rcds[j]);
					}
					else {
						ldd2[ldd2_n].Add(bit.blk.rcds[j]);
						if(ldd2[ldd2_n].n == 100){
							++ldd2_n;
							if(ldd2_n == ldd2_m){
								twk1_block_t* temp = ldd2;
								ldd2 = new twk1_block_t[ldd2_m*2];
								for(int z = 0; z < ldd2_m; ++z)
									ldd2[z] = std::move(temp[z]);
								delete[] temp;
								ldd2_m *= 2;
							}
						}
					}
				} else {
					bool found = false;
					for(int p = 0; p < mivals.size(); ++p){
						if(mivals[p].value == 0){
							if(ldd2[0].n) assert(ldd2[0].rcds[0].pos != bit.blk.rcds[j].pos);
							ldd2[0].Add(bit.blk.rcds[j]);
							found = true;
						}
					}
					if(found == false){
						std::cerr << utility::timestamp("ERROR") << "Corrupted intervals! (Too many matches)" << std::endl;
						return false;
					}
				}
			}
		}
	}

	if(ldd2[0].n == 0){
		std::cerr << "no data found for reference" << std::endl;
		return false;
	}

	if(ldd2_n == 1){
		std::cerr << "no surrounding variants" << std::endl;
		return false;
	}

	this->n_blks = ldd2_n;
	ldd = new twk1_ldd_blk[n_blks];
	for(int i = 0; i < n_blks; ++i){
		//std::cerr << "ldd2 size=" << ldd2[i].n << "/" << ldd2[i].m << " " << ldd2[i].minpos << "-" << ldd2[i].maxpos << std::endl;
		ldd[i].SetOwn(ldd2[i], reader.hdr.GetNumberSamples());
		ldd[i].Inflate(reader.hdr.GetNumberSamples(),settings.ldd_load_type, true);
		//std::cerr << ldd[i].n_rec << std::endl;
		//std::cerr << "first=" << ldd2[i].rcds[0].pos << std::endl;
	}

	return true;
}

bool twk_ld::twk_ld_impl::BuildIntervals(twk_reader& reader, const twk_ld_settings& settings, twk1_blk_iterator& bit)
{
	if(settings.ival_strings.size() == 0){ // if have interval strings
		return false;
	}

	this->intervals.ivecs.resize(reader.hdr.GetNumberContigs());
	// Parse the literal interval strings into actual interval tuples.
	if(this->ParseIntervalStrings(reader, settings) == false)
		return false;

	// Construct interval trees and find overlaps.
	if(this->intervals.Build(reader.hdr.GetNumberContigs(), reader.index) == false){
		return false;
	}

	// Number of blocks available to iterate over.
	this->n_blks = this->intervals.overlap_blocks.size();

	return true;
}

bool twk_ld::twk_ld_impl::LoadTargetBlocks(twk_reader& reader,
		twk1_blk_iterator& bit,
		const twk_ld_balancer& balancer,
		const twk_ld_settings& settings)
{
	if(intervals.overlap_blocks.size() == 0){ // if have interval strings
		return false;
	}

	// src1 and src2 blocks are the same: e.g. (5,5) or (10,10).
	if(balancer.diag){
		//std::cerr << "is diag" << std::endl;
		//std::cerr << "load range=" << balancer.toR - balancer.fromR << std::endl;

		n_blks = balancer.toR - balancer.fromR;
		m_blks = balancer.toR - balancer.fromR;
		ldd    = new twk1_ldd_blk[m_blks];
		ldd2   = new twk1_block_t[m_blks];
	}
	// Computation is square: e.g. (5,6) or (10,21).
	else {
		//std::cerr << "is square" << std::endl;
		//std::cerr << "load range=" << (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR) << std::endl;

		n_blks = (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR);
		m_blks = (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR);
		ldd    = new twk1_ldd_blk[m_blks];
		ldd2   = new twk1_block_t[m_blks];
	}

	if(settings.low_memory)
		std::cerr << utility::timestamp("LOG") << "Running in restriced memory mode..." << std::endl;
	else
		std::cerr << utility::timestamp("LOG") << "Running in standard mode. Pre-computing data..." << std::endl;

#if SIMD_AVAILABLE == 1
		std::cerr << utility::timestamp("LOG","SIMD") << "Vectorized instructions available: " << TWK_SIMD_MAPPING[SIMD_VERSION] << "..." << std::endl;
#else
		std::cerr << utility::timestamp("LOG","SIMD") << "No vectorized instructions available..." << std::endl;
#endif
		std::cerr << utility::timestamp("LOG") << "Constructing list, vector, RLE... ";

		Timer timer; timer.Start();
		if(balancer.diag){
			bit.stream->seekg(intervals.overlap_blocks[balancer.fromL]->foff);
			for(int i = 0; i < (balancer.toL - balancer.fromL); ++i){
				if(bit.NextBlock() == false){
					std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
					return false;
				}

				ldd2[i] = std::move(bit.blk);
				ldd[i].SetOwn(ldd2[i], reader.hdr.GetNumberSamples());
				ldd[i].Inflate(reader.hdr.GetNumberSamples(), settings.ldd_load_type, true);
			}
		} else {
			uint32_t offset = 0;
			bit.stream->seekg(intervals.overlap_blocks[balancer.fromL]->foff);
			for(int i = 0; i < (balancer.toL - balancer.fromL); ++i){
				if(bit.NextBlock() == false){
					std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
					return false;
				}

				ldd2[offset] = std::move(bit.blk);
				ldd[offset].SetOwn(ldd2[offset], reader.hdr.GetNumberSamples());
				ldd[offset].Inflate(reader.hdr.GetNumberSamples(),settings.ldd_load_type, true);
				++offset;
			}

			bit.stream->seekg(intervals.overlap_blocks[balancer.fromR]->foff);
			for(int i = 0; i < (balancer.toR - balancer.fromR); ++i){
				if(bit.NextBlock() == false){
					std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
					return false;
				}

				ldd2[offset] = std::move(bit.blk);
				ldd[offset].SetOwn(ldd2[offset], reader.hdr.GetNumberSamples());
				ldd[offset].Inflate(reader.hdr.GetNumberSamples(),settings.ldd_load_type, true);
				++offset;
			}
		}

		std::cerr << "Done! " << timer.ElapsedString() << std::endl;
		//std::cerr << "ldd2=" << n_blks << "/" << m_blks << std::endl;
		return(true);

	return true;
}

bool twk_ld::twk_ld_impl::LoadAllBlocks(twk_reader& reader,
		twk1_blk_iterator& bit,
		const twk_ld_balancer& balancer,
		const twk_ld_settings& settings)
{
	if(reader.index.n == 0)
		return false;

	// Delete old
	delete[] ldd; delete[] ldd2;
	ldd = nullptr; ldd2 = nullptr;

	if(balancer.diag){
		//std::cerr << "is diag" << std::endl;
		//std::cerr << "load range=" << balancer.toR - balancer.fromR << std::endl;

		n_blks = balancer.toR - balancer.fromR;
		m_blks = balancer.toR - balancer.fromR;
		std::cerr << utility::timestamp("LOG") << "Allocating " << utility::ToPrettyString(m_blks) << " blocks..." << std::endl;

		ldd  = new twk1_ldd_blk[m_blks];
		ldd2 = new twk1_block_t[m_blks];
	} else {
		//std::cerr << "is square" << std::endl;
		//std::cerr << "load range=" << (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR) << std::endl;

		n_blks = (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR);
		m_blks = (balancer.toL - balancer.fromL) + (balancer.toR - balancer.fromR);
		std::cerr << utility::timestamp("LOG") << "Allocating " << utility::ToPrettyString(m_blks) << " blocks..." << std::endl;

		ldd  = new twk1_ldd_blk[m_blks];
		ldd2 = new twk1_block_t[m_blks];
	}

	if(settings.low_memory)
		std::cerr << utility::timestamp("LOG") << "Running in restriced memory mode..." << std::endl;
	else
		std::cerr << utility::timestamp("LOG") << "Running in standard mode. Pre-computing data..." << std::endl;

#if SIMD_AVAILABLE == 1
	std::cerr << utility::timestamp("LOG","SIMD") << "Vectorized instructions available: " << TWK_SIMD_MAPPING[SIMD_VERSION] << "..." << std::endl;
#else
	std::cerr << utility::timestamp("LOG","SIMD") << "No vectorized instructions available..." << std::endl;
#endif

	// Distribute unpacking across multiple thread slaves.
	// This has insignificant performance impact on small files and/or
	// files with small number of samples (n<10,000). It has considerable
	// impact of large files when hundreds to thousands of samples and/or
	// a large number of variant blocks.
	const uint32_t rangeL = balancer.toL - balancer.fromL;
	const uint32_t rangeR = balancer.toR - balancer.fromR;
	uint32_t unpack_threads = settings.n_threads;
	uint32_t ppthreadL = std::ceil((float)(balancer.toL - balancer.fromL) / unpack_threads);
	// Assert blocks are balanced.
	assert(balancer.toL - balancer.fromL == balancer.toR - balancer.fromR);

	if(ppthreadL == 0){
		unpack_threads = rangeL;
		//std::cerr << "changing unpack threads=" << settings.n_threads << "->" << unpack_threads << std::endl;
	} else {
		if(ppthreadL * unpack_threads > rangeL){
			unpack_threads = rangeL / ppthreadL;
			//std::cerr << "changing unpack threads=" << settings.n_threads << "->" << unpack_threads << std::endl;
		}
	}

	//std::cerr << "balance=" << (balancer.toL - balancer.fromL) << " and " << (balancer.toL - balancer.fromL) / unpack_threads << " -> " << ppthreadL << std::endl;
	std::cerr << utility::timestamp("LOG") << "Constructing list, vector, RLE..." << std::endl;
	std::cerr << utility::timestamp("LOG","THREAD") << "Unpacking using " << unpack_threads << " threads: ";

	Timer timer; timer.Start();
	twk_ld_unpacker* slaves = new twk_ld_unpacker[unpack_threads];
	for(uint32_t i = 0; i < unpack_threads; ++i){
		slaves[i].ldd  = ldd;
		slaves[i].ldd2 = ldd2;
		slaves[i].rdr  = &reader;
		slaves[i].resize = true;
		slaves[i].fL = balancer.fromL + ppthreadL*i; // from-left
		slaves[i].tL = i+1 == unpack_threads ? balancer.fromL+rangeL : balancer.fromL+(ppthreadL*(i+1)); // to-left
		slaves[i].fR = balancer.fromR + ppthreadL*i; // from-right
		slaves[i].tR = i+1 == unpack_threads ? balancer.fromR+rangeR : balancer.fromR+(ppthreadL*(i+1)); // to-right
		slaves[i].loff   = balancer.toL - balancer.fromL; // offset from left
		slaves[i].lshift = slaves[i].fL - balancer.fromL; // offset from right
		slaves[i].roff   = slaves[i].fR - balancer.fromR;
		std::cerr << '.'; std::cerr.flush();
		//std::cerr << "range=" << slaves[i].fL << "->" << slaves[i].tL << " and " << slaves[i].fR << "->" << slaves[i].tR << std::endl;
	}

	for(uint32_t i = 0; i < unpack_threads; ++i) slaves[i].Start(balancer.diag, settings.ldd_load_type, settings.in);
	for(uint32_t i = 0; i < unpack_threads; ++i) slaves[i].thread->join();
	delete[] slaves;

	std::cerr << " Done! " << timer.ElapsedString() << std::endl;
	return(true);
}

bool twk_ld::Compute(const twk_ld_settings& settings){
	this->settings = settings;
	return(Compute());
}

bool twk_ld::ComputeSingle(const twk_ld_settings& settings, bool verbose, bool progress){
	this->settings = settings;
	return(ComputeSingle(verbose, progress));
}

bool twk_ld::Compute(){
	//return(ComputePerformance());

	if(settings.in.size() == 0){
		std::cerr << utility::timestamp("ERROR") << "No file-name provided..." << std::endl;
		return false;
	}

	if(settings.window && settings.n_chunks != 1){
		std::cerr << utility::timestamp("ERROR") << "Cannot use chunking in window mode!" << std::endl;
		return false;
	}

	std::cerr << utility::timestamp("LOG","READER") << "Opening " << settings.in << "..." << std::endl;

	twk_reader reader;
	if(reader.Open(settings.in) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to open file: " << settings.in << "..." << std::endl;
		return false;
	}

	std::cerr << utility::timestamp("LOG") << "Samples: " << utility::ToPrettyString(reader.hdr.GetNumberSamples()) << "..." << std::endl;

	twk1_blk_iterator bit;
	bit.stream = reader.stream;

	Timer timer;

	/**<
	 * Trigger the appropriate flag(s) for constructing bitvector, bitmap, or
	 * list structures. These additional structs are what we actually compute
	 * LD from.
	 */
	if(settings.low_memory && settings.force_phased && settings.bitmaps){
		settings.ldd_load_type = TWK_LDD_BITMAP;
	} else if(settings.low_memory && settings.force_phased){
		settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
	} else if(settings.force_phased){
		settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
	} else if(settings.forced_unphased){
		settings.ldd_load_type = TWK_LDD_VEC;
	} else {
		settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
	}

	// Construct interval trees for the input interval strings, if any.
	if(settings.ival_strings.size() == 0) mImpl->n_blks = reader.index.n;
	else {
		if(this->mImpl->BuildIntervals(reader, settings, bit) == false)
			return false;
	}

	if(settings.window) settings.c_chunk = 0;

	twk_ld_balancer balancer;
	if(balancer.Build(mImpl->n_blks, settings.n_chunks, settings.c_chunk) == false){
		return false;
	}
	std::cerr << utility::timestamp("LOG","BALANCING") << "Using ranges [" << balancer.fromL << "-" << balancer.toL << "," << balancer.fromR << "-" << balancer.toR << "] in " << (settings.window ? "window mode" : "square mode") <<"..." << std::endl;

	// Load and construct data blocks.
	if(this->mImpl->LoadBlocks(reader, bit, balancer, settings) == false){
		return false;
	}

	if(mImpl->n_blks == 0){
		std::cerr << utility::timestamp("ERROR") << "No valid data available..." << std::endl;
		return true;
	}

	// Compute workload for progress monitoring.
	uint32_t n_variants = 0;
	uint64_t n_comparisons = 0;
	if(balancer.diag){
		for(int i = balancer.fromL; i < balancer.toL; ++i) n_variants += reader.index.ent[i].n;
		n_comparisons = ((uint64_t)n_variants * n_variants - n_variants) / 2;
		std::cerr << utility::timestamp("LOG") << utility::ToPrettyString(n_variants) << " variants from " << utility::ToPrettyString(balancer.n_m) << " blocks..." << std::endl;
	} else {
		uint32_t n_variants_left = 0, n_variants_right = 0;
		for(int i = balancer.fromL; i < balancer.toL; ++i){
			n_variants += reader.index.ent[i].n;
			n_variants_left += reader.index.ent[i].n;
		}
		for(int i = balancer.fromR; i < balancer.toR; ++i){
			n_variants += reader.index.ent[i].n;
			n_variants_left += reader.index.ent[i].n;
		}
		n_comparisons = n_variants_left * n_variants_right;
		std::cerr << utility::timestamp("LOG") << utility::ToPrettyString(n_variants) << " variants (" << utility::ToPrettyString(n_variants_left) << "," << utility::ToPrettyString(n_variants_right) << ") from " << utility::ToPrettyString(balancer.n_m) << " blocks..." << std::endl;
	}
	std::cerr << utility::timestamp("LOG","PARAMS") << settings.GetString() << std::endl;
	std::cerr << utility::timestamp("LOG") << "Performing: " << utility::ToPrettyString(n_comparisons) << " variant comparisons..." << std::endl;

	twk_ld_dynamic_balancer ticker;
	ticker = balancer;
	ticker.SetWindow(settings.window, settings.l_window);
	ticker.ldd  = mImpl->ldd;

	twk_ld_progress progress;
	progress.n_s = reader.hdr.GetNumberSamples();
	if(settings.window == false){
		progress.n_cmps = n_comparisons;
	}
	twk_ld_slave* slaves = new twk_ld_slave[settings.n_threads];
	std::vector<std::thread*> threads(settings.n_threads);

	// Start writing file.
	twk_two_writer_t* writer = nullptr;
	if(settings.out.size() == 0 || (settings.out.size() == 1 && settings.out[0] == '-')){
		std::cerr << utility::timestamp("LOG","WRITER") << "Writing to " << "stdout..." << std::endl;
		writer = new twk_two_writer_t;
	} else {
		std::string base_path = twk_writer_t::GetBasePath(settings.out);
		std::string base_name = twk_writer_t::GetBaseName(settings.out);
		std::string extension = twk_writer_t::GetExtension(settings.out);
		if(extension.length() == 3){
			if(strncasecmp(&extension[0], "two", 3) != 0){
				settings.out =  (base_path.size() ? base_path + "/" : "") + base_name + ".two";
			}
		} else {
			 settings.out = (base_path.size() ? base_path + "/" : "") + base_name + ".two";
		}

		std::cerr << utility::timestamp("LOG","WRITER") << "Opening " << settings.out << "..." << std::endl;
		writer = new twk_two_writer_t;
		if(writer->Open(settings.out) == false){
			std::cerr << utility::timestamp("ERROR", "WRITER") << "Failed to open file: " << settings.out << "..." << std::endl;
			delete writer;
			return false;
		}
	}

	// Append literal string.
	std::string calc_string = "\n##tomahawk_calcVersion=" + std::string(VERSION) + "\n";
	calc_string += "##tomahawk_calcCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + utility::datetime() + "\n";
	reader.hdr.literals_ += calc_string;

	if(writer->WriteHeaderBinary(reader) == false){
		std::cerr << utility::timestamp("ERROR","WRITER") << "Failed to write header!" << std::endl;
		return false;
	}
	// end start write

	// New index
	IndexOutput index(reader.hdr.GetNumberContigs());

	timer.Start();
	std::cerr << utility::timestamp("LOG","THREAD") << "Spawning " << settings.n_threads << " threads: ";
	for(int i = 0; i < settings.n_threads; ++i){
		slaves[i].ldd    = mImpl->ldd;
		slaves[i].n_s    = reader.hdr.GetNumberSamples();
		slaves[i].ticker = &ticker;
		slaves[i].engine.SetSamples(reader.hdr.GetNumberSamples());
		slaves[i].engine.SetBlocksize(settings.b_size);
		slaves[i].engine.progress = &progress;
		slaves[i].engine.writer   = writer;
		slaves[i].engine.index    = &index;
		slaves[i].engine.settings = settings;
		slaves[i].progress = &progress;
		slaves[i].settings = &settings;
		threads[i] = slaves[i].Start();
		std::cerr << ".";
	}
	std::cerr << std::endl;

	progress.Start();

	for(int i = 0; i < settings.n_threads; ++i) threads[i]->join();
	for(int i = 0; i < settings.n_threads; ++i) slaves[i].engine.CompressBlock();
	progress.is_ticking = false;
	progress.PrintFinal();
	writer->stream.flush();

	/*
	std::cerr << utility::timestamp("LOG","THREAD") << "Thread\tOutput\tTWK-LIST\tTWK-BVP-BM\tTWK-BVP\tTWK-BVP-NM\tTWK-BVU\tTWK-BVU-NM\tTWK-RLEP\tTWK-RLEU\n";
	for(int i = 0; i < settings.n_threads; ++i){
		std::cerr << i << "\t" << utility::ToPrettyString(slaves[i].engine.n_out);
		for(int j = 0; j < 8; ++j){
			std::cerr << "\t" << utility::ToPrettyString(slaves[i].engine.n_method[j]);
		}
		std::cerr << std::endl;
	}
	*/

	//std::cerr << "performed=" << ticker.n_perf << std::endl;
	if(writer->WriteFinal(index) == false){
		std::cerr << utility::timestamp("ERROR","WRITER") << "Failed to write final block!" << std::endl;
		return false;
	}

	delete[] slaves; delete writer;
	std::cerr << utility::timestamp("LOG","PROGRESS") << "All done..." << timer.ElapsedString() << "!" << std::endl;

	return true;
}

bool twk_ld::ComputeSingle(bool verbose, bool progress){
	if(settings.in.size() == 0){
		std::cerr << utility::timestamp("ERROR") << "No file-name provided..." << std::endl;
		return false;
	}

	if(settings.n_chunks != 1){
		std::cerr << utility::timestamp("ERROR") << "Cannot use chunking in single mode!" << std::endl;
		return false;
	}

	if(settings.window){
		std::cerr << utility::timestamp("ERROR") << "Cannot use window in single mode!" << std::endl;
		return false;
	}

	if(settings.ival_strings.size() == 0){
		std::cerr << utility::timestamp("ERROR") << "An interval has to be provided in single mode!" << std::endl;
		return false;
	}

	if(settings.ival_strings.size() != 1){
		std::cerr << utility::timestamp("ERROR") << "Only a single interval can be provided in single mode!" << std::endl;
		return false;
	}

	if(verbose) std::cerr << utility::timestamp("LOG","READER") << "Opening " << settings.in << "..." << std::endl;

	twk_reader reader;
	if(reader.Open(settings.in) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to open file: " << settings.in << "..." << std::endl;
		return false;
	}

	if(verbose) std::cerr << utility::timestamp("LOG") << "Samples: " << utility::ToPrettyString(reader.hdr.GetNumberSamples()) << "..." << std::endl;

	twk1_blk_iterator bit;
	bit.stream = reader.stream;

	Timer timer;

	/**<
	 * Trigger the appropriate flag(s) for constructing bitvector, bitmap, or
	 * list structures. These additional structs are what we actually compute
	 * LD from.
	 */
	if(settings.low_memory && settings.force_phased && settings.bitmaps){
		settings.ldd_load_type = TWK_LDD_BITMAP;
	} else if(settings.low_memory && settings.force_phased){
		settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
	} else if(settings.force_phased){
		settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
	} else if(settings.forced_unphased){
		settings.ldd_load_type = TWK_LDD_VEC;
	} else {
		settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
	}

	// Construct interval trees for the input interval string.
	if(this->mImpl->LoadTargetSingle(reader, bit, settings, settings.ldd_load_type) == false)
		return false;

	twk_ld_balancer balancer;
	if(balancer.BuildSingleSite(mImpl->n_blks, settings.n_chunks, settings.c_chunk) == false)
		return false;

	if(verbose) std::cerr << utility::timestamp("LOG","BALANCING") << "Using ranges [" << balancer.fromL << "-" << balancer.toL << "," << balancer.fromR << "-" << balancer.toR << "] in " << (settings.window ? "window mode" : "square mode") <<"..." << std::endl;

	if(mImpl->n_blks == 0){
		std::cerr << utility::timestamp("ERROR") << "No valid data available..." << std::endl;
		return true;
	}

	// Compute workload for progress monitoring.
	// Todo
	uint32_t n_variants = 0;
	uint64_t n_comparisons = 0;
	if(balancer.diag){
		for(int i = balancer.fromL; i < balancer.toL; ++i) n_variants += reader.index.ent[i].n;
		n_comparisons = ((uint64_t)n_variants * n_variants - n_variants) / 2;
		if(verbose) std::cerr << utility::timestamp("LOG") << utility::ToPrettyString(n_variants) << " variants from " << utility::ToPrettyString(balancer.n_m) << " blocks..." << std::endl;
	} else {
		uint32_t n_variants_left = 0, n_variants_right = 0;
		for(int i = balancer.fromL; i < balancer.toL; ++i){
			n_variants += reader.index.ent[i].n;
			n_variants_left += reader.index.ent[i].n;
		}
		for(int i = balancer.fromR; i < balancer.toR; ++i){
			n_variants += reader.index.ent[i].n;
			n_variants_left += reader.index.ent[i].n;
		}
		n_comparisons = n_variants_left * n_variants_right;
		if(verbose) std::cerr << utility::timestamp("LOG") << utility::ToPrettyString(n_variants) << " variants (" << utility::ToPrettyString(n_variants_left) << "," << utility::ToPrettyString(n_variants_right) << ") from " << utility::ToPrettyString(balancer.n_m) << " blocks..." << std::endl;
	}

	if(verbose){
		std::cerr << utility::timestamp("LOG","PARAMS") << settings.GetString() << std::endl;
		std::cerr << utility::timestamp("LOG") << "Performing: " << utility::ToPrettyString(n_comparisons) << " variant comparisons..." << std::endl;
	}

	twk_ld_dynamic_balancer ticker;
	ticker = balancer;
	ticker.SetWindow(settings.window, settings.l_window);
	ticker.ldd = mImpl->ldd;

	twk_ld_progress progression;
	progression.n_s = reader.hdr.GetNumberSamples();
	if(settings.window == false){
		progression.n_cmps = n_comparisons;
	}
	twk_ld_slave* slaves = new twk_ld_slave[settings.n_threads];
	std::vector<std::thread*> threads(settings.n_threads);

	// Start writing file.
	twk_two_writer_t* writer = nullptr;
	if(settings.out.size() == 0 || (settings.out.size() == 1 && settings.out[0] == '-')){
		if(verbose) std::cerr << utility::timestamp("LOG","WRITER") << "Writing to " << "stdout..." << std::endl;
		writer = new twk_two_writer_t;
	} else {
		std::string base_path = twk_writer_t::GetBasePath(settings.out);
		std::string base_name = twk_writer_t::GetBaseName(settings.out);
		std::string extension = twk_writer_t::GetExtension(settings.out);
		if(extension.length() == 3){
			if(strncasecmp(&extension[0], "two", 3) != 0){
				settings.out =  (base_path.size() ? base_path + "/" : "") + base_name + ".two";
			}
		} else {
			 settings.out = (base_path.size() ? base_path + "/" : "") + base_name + ".two";
		}

		if(verbose) std::cerr << utility::timestamp("LOG","WRITER") << "Opening " << settings.out << "..." << std::endl;
		writer = new twk_two_writer_t;
		if(writer->Open(settings.out) == false){
			std::cerr << utility::timestamp("ERROR", "WRITER") << "Failed to open file: " << settings.out << "..." << std::endl;
			delete writer;
			return false;
		}
	}

	// Append literal string.
	std::string calc_string = "\n##tomahawk_calcVersion=" + std::string(VERSION) + "\n";
	calc_string += "##tomahawk_calcCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + utility::datetime() + "\n";
	reader.hdr.literals_ += calc_string;

	if(writer->WriteHeaderBinary(reader) == false){
		std::cerr << utility::timestamp("ERROR","WRITER") << "Failed to write header!" << std::endl;
		return false;
	}
	// end start write

	// New index
	IndexOutput index(reader.hdr.GetNumberContigs());

	timer.Start();
	if(verbose) std::cerr << utility::timestamp("LOG","THREAD") << "Spawning " << settings.n_threads << " threads: ";
	for(int i = 0; i < settings.n_threads; ++i){
		slaves[i].ldd    = mImpl->ldd;
		slaves[i].n_s    = reader.hdr.GetNumberSamples();
		slaves[i].ticker = &ticker;
		slaves[i].engine.SetSamples(reader.hdr.GetNumberSamples());
		slaves[i].engine.SetBlocksize(settings.b_size);
		slaves[i].engine.progress = &progression;
		slaves[i].engine.writer   = writer;
		slaves[i].engine.index    = &index;
		slaves[i].engine.settings = settings;
		slaves[i].progress = &progression;
		slaves[i].settings = &settings;
		threads[i] = slaves[i].Start();
		if(verbose) std::cerr << ".";
	}
	if(verbose) std::cerr << std::endl;

	if(progress) progression.Start();

	for(int i = 0; i < settings.n_threads; ++i) threads[i]->join();
	for(int i = 0; i < settings.n_threads; ++i) slaves[i].engine.CompressBlock();
	progression.is_ticking = false;
	if(progress) progression.PrintFinal();
	writer->stream.flush();

	/*
	if(verbose){
		std::cerr << utility::timestamp("LOG","THREAD") << "Thread\tOutput\tTWK-LIST\tTWK-BVP-BM\tTWK-BVP\tTWK-BVP-NM\tTWK-BVU\tTWK-BVU-NM\tTWK-RLEP\tTWK-RLEU\n";
		for(int i = 0; i < settings.n_threads; ++i){
			std::cerr << i << "\t" << utility::ToPrettyString(slaves[i].engine.n_out);
			for(int j = 0; j < 8; ++j){
				std::cerr << "\t" << utility::ToPrettyString(slaves[i].engine.n_method[j]);
			}
			std::cerr << std::endl;
		}
	}
	*/

	//std::cerr << "performed=" << ticker.n_perf << std::endl;
	if(writer->WriteFinal(index) == false){
		std::cerr << utility::timestamp("ERROR","WRITER") << "Failed to write final block!" << std::endl;
		return false;
	}

	delete[] slaves; delete writer;
	if(verbose) std::cerr << utility::timestamp("LOG","PROGRESS") << "All done..." << timer.ElapsedString() << "!" << std::endl;

	return true;
}

bool twk_ld::ComputePerformance(){
	if(TWK_SLAVE_DEBUG_MODE != 1){
		std::cerr << utility::timestamp("ERROR","COMPILATION") << "Cannot run performance mode without compiling SLAVE_DEBUG_MODE set to 1!" << std::endl;
		return false;
	}

	if(settings.in.size() == 0){
		std::cerr << utility::timestamp("ERROR") << "No file-name provided..." << std::endl;
		return false;
	}

	if(settings.window && settings.n_chunks != 1){
		std::cerr << utility::timestamp("ERROR") << "Cannot use chunking in window mode!" << std::endl;
		return false;
	}

	std::cerr << utility::timestamp("LOG","READER") << "Opening " << settings.in << "..." << std::endl;

	twk_reader reader;
	if(reader.Open(settings.in) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to open file: " << settings.in << "..." << std::endl;
		return false;
	}

	std::cerr << utility::timestamp("LOG") << "Samples: " << utility::ToPrettyString(reader.hdr.GetNumberSamples()) << "..." << std::endl;

	twk1_blk_iterator bit;
	bit.stream = reader.stream;

	Timer timer;

	//settings.ldd_load_type = TWK_LDD_ALL;
	if(settings.low_memory && settings.force_phased && settings.bitmaps){
		settings.ldd_load_type = TWK_LDD_BITMAP;
	} else if(settings.low_memory && settings.force_phased){
		settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
	} else if(settings.force_phased){
		settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
	} else if(settings.forced_unphased){
		settings.ldd_load_type = TWK_LDD_VEC;
	} else {
		settings.ldd_load_type = TWK_LDD_VEC | TWK_LDD_LIST;
	}

	if(settings.ival_strings.size() == 0) mImpl->n_blks = reader.index.n;
	else {
		if(this->mImpl->BuildIntervals(reader, settings, bit) == false)
			return false;
	}

	if(settings.window) settings.c_chunk = 0;

	twk_ld_balancer balancer;
	if(balancer.Build(mImpl->n_blks, settings.n_chunks, settings.c_chunk) == false){
		return false;
	}

	std::cerr << utility::timestamp("LOG","BALANCING") << "Using ranges [" << balancer.fromL << "-" << balancer.toL << "," << balancer.fromR << "-" << balancer.toR << "] in " << (settings.window ? "window mode" : "square mode") <<"..." << std::endl;

	if(mImpl->n_blks == 0){
		std::cerr << utility::timestamp("ERROR") << "No valid data available..." << std::endl;
		return true;
	}

	uint32_t n_variants = 0;
	if(balancer.diag){
		for(int i = balancer.fromL; i < balancer.toL; ++i) n_variants += reader.index.ent[i].n;
	} else {
		for(int i = balancer.fromL; i < balancer.toL; ++i) n_variants += reader.index.ent[i].n;
		for(int i = balancer.fromR; i < balancer.toR; ++i) n_variants += reader.index.ent[i].n;
	}

	std::cerr << utility::timestamp("LOG") << utility::ToPrettyString(n_variants) << " variants from " << utility::ToPrettyString(balancer.n_m) << " blocks..." << std::endl;
	std::cerr << utility::timestamp("LOG","PARAMS") << settings.GetString() << std::endl;
	std::cerr << utility::timestamp("LOG") << "Performing: " << utility::ToPrettyString(((uint64_t)n_variants * n_variants - n_variants) / 2) << " variant comparisons..." << std::endl;

	twk_ld_progress progress;
	progress.n_s = reader.hdr.GetNumberSamples();
	if(settings.window == false){
		if(balancer.diag)
			progress.n_cmps = ((uint64_t)n_variants * n_variants - n_variants) / 2;
		else {

		}
	}

	// Append literal string.
	std::string calc_string = "\n##tomahawk_calcVersion=" + std::string(VERSION) + "\n";
	calc_string += "##tomahawk_calcCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + utility::datetime() + "\n";
	reader.hdr.literals_ += calc_string;

	// New index
	IndexOutput index(reader.hdr.GetNumberContigs());

	twk_ld_engine::ep f; // function pointer array
	f[0] = &twk_ld_engine::PhasedVectorized;
	f[1] = &twk_ld_engine::PhasedVectorizedNoMissing;
	f[2] = &twk_ld_engine::UnphasedVectorized;
	f[3] = &twk_ld_engine::UnphasedVectorizedNoMissing;
	f[4] = &twk_ld_engine::PhasedRunlength;
	f[5] = &twk_ld_engine::PhasedList;
	f[6] = &twk_ld_engine::UnphasedRunlength;
	f[7] = &twk_ld_engine::PhasedBitmap;
	//f[8] = &twk_ld_engine::PhasedListSpecial;
	f[8] = &twk_ld_engine::PhasedListVector;

	const std::vector<std::string> method_names = {"PhasedVectorized",
			"PhasedVectorizedNoMissing",
			"UnphasedVectorized",
			"UnphasedVectorizedNoMissing",
			"PhasedRunlength",
			"PhasedList",
			"UnphasedRunlength",
			"PhasedBitmap",
			"PhasedListSpecial",
			"PhasedListBV"};

	uint8_t unpack_list[9];
	memset(unpack_list, TWK_LDD_VEC, 9);
	unpack_list[4] = TWK_LDD_NONE;
	unpack_list[5] = TWK_LDD_LIST;
	unpack_list[6] = TWK_LDD_NONE;
	unpack_list[7] = TWK_LDD_BITMAP;
	unpack_list[8] = TWK_LDD_VEC | TWK_LDD_LIST;
	//unpack_list[9] = TWK_LDD_VEC | TWK_LDD_LIST;

	twk_ld_perf perfs[10];
	uint32_t perf_size = reader.hdr.GetNumberSamples()*4 + 2;
	for(int i = 0; i < 9; ++i){
		perfs[i].cycles = new uint64_t[perf_size];
		perfs[i].freq   = new uint64_t[perf_size];
		memset(perfs[i].cycles, 0, perf_size*sizeof(uint64_t));
		memset(perfs[i].freq,   0, perf_size*sizeof(uint64_t));
	}

	//prepare_shuffling_dictionary();
	for(int method = 8; method < 9; ++method){
		if(method == 2 || method == 3 || method == 4 || method == 6) continue;
		std::cerr << utility::timestamp("LOG","PROGRESS") << "Starting method: " << method_names[method] << "..." << std::endl;
		settings.ldd_load_type = unpack_list[method];
		if(this->mImpl->LoadBlocks(reader, bit, balancer, settings) == false){
			return false;
		}

		twk_ld_dynamic_balancer ticker;
		ticker = balancer;
		ticker.SetWindow(settings.window, settings.l_window);
		ticker.ldd  = mImpl->ldd;

		timer.Start();
		twk_ld_slave s;
		s.ldd    = mImpl->ldd;
		s.n_s    = reader.hdr.GetNumberSamples();
		s.ticker = &ticker;
		s.engine.SetSamples(reader.hdr.GetNumberSamples());
		s.engine.SetBlocksize(settings.b_size);
		s.engine.progress = &progress;
		//s.engine.writer   = writer;
		s.engine.index    = &index;
		s.engine.settings = settings;
		s.progress = &progress;
		s.settings = &settings;

		s.CalculatePerformance(f[method], &perfs[method]);

		std::cerr << utility::timestamp("LOG","PROGRESS") << "Finished method: " << method_names[method] << ". Elapsed time=" << timer.ElapsedString() << std::endl;
		//std::cerr << "m=" << method << " " << utility::ToPrettyString(1) << " " << utility::ToPrettyString((uint64_t)((float)1/timer.Elapsed().count())) << "/s " << utility::ToPrettyString((uint64_t)((float)1*reader.hdr.GetNumberSamples()/timer.Elapsed().count())) << "/gt/s " << timer.ElapsedString() << " " << utility::ToPrettyString(1) << std::endl;
	}

	for(int i = 0; i < 9; ++i){
		for(int j = 0; j < perf_size; ++j){
			// Print if non-zero
			if(perfs[i].freq[j] != 0)
				std::cout << i << "\t" << j << "\t" << (double)perfs[i].cycles[j]/perfs[i].freq[j] << "\t" << perfs[i].freq[j] << "\t" << perfs[i].cycles[j] << '\n';
		}
		delete[] perfs[i].cycles;
	}

	return true;
}

}
