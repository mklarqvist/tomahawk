#ifndef LIB_TWO_SORTER_STRUCTS_H_
#define LIB_TWO_SORTER_STRUCTS_H_

#include "sort_progress.h"
#include "zstd_codec.h"
#include "two_reader.h"

namespace tomahawk {

struct sort_helper {
	sort_helper() : rid(0), minP(0), maxP(0), foff(0), fend(0), n(0), nc(0){}

	uint32_t rid, minP, maxP;
	uint64_t foff, fend, n, nc;
};

struct twk_sort_slave {
	struct run_intervals {
		uint32_t ref_rid, n_run, minp, maxp;
	};

	std::thread* Start(IndexOutput& rdr);

	bool Sort();

public:
	float m_limit;
	uint32_t n; // limit in gb, number of entries that corresponds to
	uint32_t f, t, bl_size, c_level; // (from,to)-tuple, flush block-size
	std::string filename, tmp_filename; // temporary filename
	std::thread* thread;
	std::ifstream stream;
	std::ofstream ostream;
	twk1_two_iterator* it;
	twk1_two_block_t* blk;
	std::vector<sort_helper> local_idx; // local offset index
	ZSTDCodec zcodec;
	twk_buffer_t obuf, obuf2;
	twk_sort_progress* progress;
	std::vector< std::vector<run_intervals> > run_ivals;
};

struct twk_two_stream_iterator {
	bool Open(const std::string file,
			const uint64_t foff,
			const uint64_t fend,
			const uint32_t n_uncompressed,
			const uint32_t n_compressed);

	bool Next(twk1_two_t& rec, const uint32_t b_read = 256000);
	bool NextBlock(const uint32_t b_read = 256000);

public:
	uint32_t n, it, n_tot, it_tot; // n_records, it_pos
	uint32_t b_left, n_unc, n_cmp, n_unc_cum;
	uint64_t off_start, off_end;
	std::ifstream stream;
	twk_buffer_t out, out2;
	ZSTDCodec zcodec;
};

struct two_queue_entry {
public:
	two_queue_entry(const twk1_two_t& data, uint32_t streamID) :
		qid(streamID), rec(data)
	{}

	inline bool operator<(const two_queue_entry& other) const {
		return(!(rec < other.rec));
	}

public:
	uint32_t qid;
	twk1_two_t rec;
};

}

#endif /* LIB_TWO_SORTER_STRUCTS_H_ */
