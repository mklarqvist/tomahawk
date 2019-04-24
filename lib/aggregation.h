#ifndef TWK_AGGREGATION_H_
#define TWK_AGGREGATION_H_

#include <thread>
#include <limits>
#include <cstdint>

#include "core.h"
#include "two_reader.h"
#include "sort_progress.h"

namespace tomahawk {

struct twk_agg_slave {
public:
	struct range_helper {
		range_helper() : set(false), min(std::numeric_limits<uint32_t>::max()), max(0){}

		bool set;
		uint32_t min, max;
	};

public:
	twk_agg_slave() : f(0), t(0), xrange(0), yrange(0), it(nullptr),
	    thread(nullptr), progress(nullptr),
	    aggregator(&twk_sstats::AddR2), reductor(&twk_sstats::GetMean)
	{}
	~twk_agg_slave(){ delete it; delete thread; }

	/**<
	 *
	 * @param rdr Reference instance of a two reader.
	 * @return    Returns a pointer to the spawned thread instance if successful or a nullptr otherwise.
	 */
	std::thread* StartFindRanges(two_reader& rdr){
		if(f > t) return nullptr;
		if(t - f == 0) return nullptr;

		if(stream.good()) stream.close();

		stream.open(filename,std::ios::binary | std::ios::in);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to open \"" << filename << "\"..." << std::endl;
			return nullptr;
		}

		stream.seekg(rdr.index.ent[f].foff);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to seek to position " << rdr.index.ent[f].foff << " in \""  << filename << "\"..." << std::endl;
			return nullptr;
		}

		delete it; it = nullptr;
		it  = new twk1_two_iterator;
		it->stream = &stream;

		// Vector of contigs.
		contig_avail.resize(rdr.hdr.GetNumberContigs());

		// New thread.
		delete thread; thread = nullptr;
		thread = new std::thread(&twk_agg_slave::FindRangesUnsorted, this);

		return(thread);
	}

	/**<
	 *
	 * @param rdr Reference instance of a two reader.
	 * @param agg Aggregation function pointer.
	 * @param red Reduction function pointer.
	 * @param x   Number of bins in X-dimension.
	 * @param y   Number of bins in Y-dimension.
	 * @param xr  Range of X.
	 * @param yr  Range of Y.
	 * @return    Returns a pointer to the spawned thread instance if successful or a nullptr otherwise.
	 */
	std::thread* StartBuildMatrix(two_reader& rdr,
			twk_sstats::aggfunc agg,
			twk_sstats::redfunc red,
			uint32_t x, uint32_t y,
			uint32_t xr, uint32_t yr)
	{
		if(f > t) return nullptr;
		if(t - f == 0) return nullptr;

		assert(agg != nullptr);
		assert(red != nullptr);
		aggregator = agg;
		reductor = red;
		xrange = xr;
		yrange = yr;
		mat = std::vector< std::vector<twk_sstats> >(x, std::vector<twk_sstats>(y));

		if(stream.good()) stream.close();

		stream.open(filename,std::ios::binary | std::ios::in);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to open \"" << filename << "\"..." << std::endl;
			return nullptr;
		}

		stream.seekg(rdr.index.ent[f].foff);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to seek to position " << rdr.index.ent[f].foff << " in \""  << filename << "\"..." << std::endl;
			return nullptr;
		}

		delete it; it = nullptr;
		it  = new twk1_two_iterator;
		it->stream = &stream;

		// New thread.
		delete thread; thread = nullptr;
		thread = new std::thread(&twk_agg_slave::BuildMatrix, this);

		return(thread);
	}

	/**<
	 * Investigate the presence of blocks mapping to unique chromosomes. This
	 * information is needed to set the (X,Y)-coordinate system for the aggregation
	 * space.
	 * @param rdr Input reference reader.
	 * @return    Returns TRUE upon success or FALSE otherwise.
	 */
	bool FindRangesUnsorted(){
		uint32_t tot = 0;
		for(int i = f; i < t; ++i){
			assert(it->NextBlock());
			// Todo: bugfix this case
			if(it->blk.n < 5) continue;

			tot += it->GetBlock().n;
			//std::cerr << "here" << std::endl;

			for(int j = 0; j < it->blk.n; ++j){
				//std::cerr << j << "/" << it->blk.n << "\t" << it->blk[j].ridA << "/" << contig_avail.size() << std::endl;
				assert(it->blk[j].ridA < contig_avail.size());
				assert(it->blk[j].ridB < contig_avail.size());

				contig_avail[it->blk[j].ridA].set = true;
				contig_avail[it->blk[j].ridB].set = true;
				contig_avail[it->blk[j].ridA].min = std::min(it->blk[j].Apos, contig_avail[it->blk[j].ridA].min);
				contig_avail[it->blk[j].ridA].max = std::max(it->blk[j].Apos, contig_avail[it->blk[j].ridA].max);
				contig_avail[it->blk[j].ridB].min = std::min(it->blk[j].Bpos, contig_avail[it->blk[j].ridB].min);
				contig_avail[it->blk[j].ridB].max = std::max(it->blk[j].Bpos, contig_avail[it->blk[j].ridB].max);
			}
			if(progress != nullptr) progress->cmps += it->GetBlock().n;
		}

		delete it; it = nullptr;
		return true;
	}

	bool BuildMatrix(){
		uint32_t tot = 0;
		for(int i = f; i < t; ++i){
			assert(it->NextBlock());
			// Todo: bugfix this case
			if(it->blk.n < 5) continue;

			tot += it->GetBlock().n;
			for(int j = 0; j < it->blk.n; ++j){
				// Invoke aggregator function.
				// Position: cumulative offset up to chromosome + left-adjusted position
				// Position: (chromosome_offset.range - chromosome_offset.max - chromoosme_offset.min) + (Apos - smallest_in_chr)
				(mat[((rid_offsets[it->blk[j].ridA].range - (rid_offsets[it->blk[j].ridA].max - rid_offsets[it->blk[j].ridA].min)) + (it->blk[j].Apos - rid_offsets[it->blk[j].ridA].min))/xrange][((rid_offsets[it->blk[j].ridB].range - (rid_offsets[it->blk[j].ridB].max - rid_offsets[it->blk[j].ridB].min)) + (it->blk[j].Bpos - rid_offsets[it->blk[j].ridB].min))/yrange].*aggregator)(&it->blk[j]);
			}
			if(progress != nullptr) progress->cmps += it->GetBlock().n;
		}

		delete it; it = nullptr;
		return true;
	}

	/**<
	 * Reduction helper for matrices. Adds the data from the matrix in the passed
	 * slave instance to the current instance.
	 * @param other Reference to other slave instance.
	 */
	void AddMatrix(const twk_agg_slave& other){
		for(int i = 0; i < mat.size(); ++i){
			for(int j = 0; j < mat[i].size(); ++j)
				mat[i][j] += other.mat[i][j];
		}
	}

	/**<
	 * Simple print subroutine for the internal matrix.
	 * @param stream     Target output stream to print to.
	 * @param min_cutoff Minimum number of elements require to output a value. If the observed frequency is smaller than this value then we emit 0.
	 */
	void PrintMatrix(std::ostream& stream, uint32_t min_cutoff = 5) const {
		for(int i = 0; i < mat.size(); ++i){
			stream << (mat[i][0].*reductor)(min_cutoff);
			for(int j = 1; j < mat[i].size(); ++j)
				stream << '\t' << (mat[i][j].*reductor)(min_cutoff);
			stream << '\n';
		}
		stream.flush();
	}

	twk1_aggregate_t& Overload(twk1_aggregate_t& agg, uint32_t min_cutoff = 5) const {
		for(int i = 0; i < mat.size(); ++i){
			agg.data[i*mat.size() + 0] = (mat[i][0].*reductor)(min_cutoff);
			for(int j = 1; j < mat[i].size(); ++j){
				agg.data[i*mat.size() + j] = (mat[i][j].*reductor)(min_cutoff);
				//stream << '\t' << (mat[i][j].*reductor)(min_cutoff);
			}
			//stream << '\n';
		}
		//stream.flush();
		return(agg);
	}

public:
	uint32_t f, t; // (from,to)-tuple
	uint32_t xrange, yrange; // (x,y)-tuple
	std::ifstream stream;
	twk1_two_iterator* it;
	std::thread* thread;
	twk_sort_progress* progress;
	twk_sstats::aggfunc aggregator; // aggregator function
	twk_sstats::redfunc reductor; // reductor function
	std::string filename; // input filename
	std::vector<range_helper> contig_avail;
	std::vector<twk1_aggregate_t::offset_tuple> rid_offsets; // mat offsets
	std::vector< std::vector<twk_sstats> > mat; // Output matrix
};

}



#endif /* LIB_AGGREGATION_H_ */
