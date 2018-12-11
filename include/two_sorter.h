#ifndef TWO_SORTER_H_
#define TWO_SORTER_H_

#include "two_reader.h"
#include "sort_progress.h"

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

	std::thread* Start(two_reader& rdr){
		if(f > t) return nullptr;
		if(t - f == 0) return nullptr;
		//delete thread;

		stream = std::ifstream(filename,std::ios::binary | std::ios::in);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to open \"" << filename << "\"..." << std::endl;
			return nullptr;
		}

		stream.seekg(rdr.index.ent[f].foff);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to seek to position " << rdr.index.ent[f].foff << " in \""  << filename << "\"..." << std::endl;
			return nullptr;
		}

		it  = new twk1_two_iterator;
		blk = new twk1_two_block_t;
		blk->resize((m_limit*1e9)/sizeof(twk1_two_t));
		it->stream = &stream;

		ostream = std::ofstream(tmp_filename, std::ios::binary | std::ios::out);
		if(ostream.good() == false){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to open temp output file \"" << tmp_filename << "\"..." << std::endl;
			return nullptr;
		}

		thread = new std::thread(&twk_sort_slave::Sort, this);
		return(thread);
	}

	bool Sort(){
		uint32_t tot = 0;
		for(int i = f; i < t; ++i){
			assert(it->NextBlock());
			tot += it->GetBlock().n;
			for(int j = 0; j < it->blk.n; ++j){
				//it->blk[j].Print(std::cerr);
				if(blk->n == blk->m){
					blk->Sort();

					//std::cerr << blk->rcds[0].ridA << " and " << blk->rcds[blk->n-1].ridA << " with n=" << blk->n << std::endl;
					//assert(blk->rcds[0].ridA == blk->rcds[blk->n-1].ridA);

					run_ivals.push_back(std::vector<run_intervals>());
					run_intervals t;
					t.ref_rid = blk->rcds[0].ridA;
					t.n_run = 1;
					t.minp = blk->rcds[0].Apos;
					t.maxp = blk->rcds[0].Apos;

					for(int p = 1; p < blk->n; ++p){
						if(blk->rcds[p].ridA != t.ref_rid){
							//std::cerr << "rid=" << t.ref_rid << ":" << t.n_run << " pos=" << t.minp << "-" << t.maxp << std::endl;
							run_ivals.back().push_back(t);
							t.n_run = 0;
							t.ref_rid = blk->rcds[p].ridA;
							t.minp = blk->rcds[p].Apos;
							t.maxp = blk->rcds[p].Apos;
						}
						++t.n_run;
						t.maxp = blk->rcds[p].Apos;
					}
					if(t.n_run){
						run_ivals.back().push_back(t);
						//std::cerr << "rid=" << t.ref_rid << ":" << t.n_run << " pos=" << t.minp << "-" << t.maxp << std::endl;
					}

					sort_helper rec;
					//rec.rid = blk->rcds[0].ridA;
					//rec.minP = blk->rcds[0].Apos;
					//rec.maxP = blk->rcds[blk->n-1].Apos;

					rec.foff = ostream.tellp();
					zcodec.InitStreamCompress(c_level);
					// Compress chunks of 10k records
					int k = 0;
					for(; k + 10000 < blk->n; k += 10000){
						for(int l = k; l < k + 10000; ++l) obuf << blk->rcds[l];
						rec.nc += zcodec.StreamCompress(obuf, obuf2, ostream, twk1_two_t::packed_size * 5000);
						progress->cmps += 10000;
						obuf.reset(); obuf2.reset();
					}

					// Compress residual records
					//std::cerr << "remainder=" << k << "/" << blk->n << std::endl;

					progress->cmps += blk->n - k;
					for(int l = k; l < blk->n; ++l) obuf << blk->rcds[l];
					rec.nc += zcodec.StreamCompress(obuf, obuf2, ostream, twk1_two_t::packed_size * 5000);
					obuf.reset(); obuf2.reset();

					zcodec.StopStreamCompress();
					//std::cerr << "after stop=" << zcodec.outbuf.pos << "/" << zcodec.outbuf.size << std::endl;
					ostream.write((const char*)zcodec.outbuf.dst, zcodec.outbuf.pos);
					rec.nc += zcodec.outbuf.pos;
					rec.n = blk->n * twk1_two_t::packed_size;
					rec.fend = ostream.tellp();
					local_idx.push_back(rec);

					//std::cerr << rec.rid << "," << rec.minP << "-" << rec.maxP << std::endl;
					blk->reset();

				}
				*blk += it->blk[j];
				//it->blk[j].Print(std::cout);
			}
		}
		// any possible remainder
		if(blk->n){
			//std::cerr << "resseting" << std::endl;
			blk->Sort();
			//obuf << blk2; // poor memory management:

			run_ivals.push_back(std::vector<run_intervals>());
			run_intervals t;
			t.ref_rid = blk->rcds[0].ridA;
			t.n_run = 1;
			t.minp = blk->rcds[0].Apos;
			t.maxp = blk->rcds[0].Apos;

			for(int p = 1; p < blk->n; ++p){
				if(blk->rcds[p].ridA != t.ref_rid){
					//std::cerr << "rid=" << t.ref_rid << ":" << t.n_run << " pos=" << t.minp << "-" << t.maxp << std::endl;
					run_ivals.back().push_back(t);
					t.n_run = 0;
					t.ref_rid = blk->rcds[p].ridA;
					t.minp = blk->rcds[p].Apos;
				}
				++t.n_run;
				t.maxp = blk->rcds[p].Apos;
			}
			if(t.n_run){
				run_ivals.back().push_back(t);
				//std::cerr << "rid=" << t.ref_rid << ":" << t.n_run << " pos=" << t.minp << "-" << t.maxp << std::endl;
			}

			sort_helper rec;
			rec.foff = ostream.tellp();
			zcodec.InitStreamCompress(c_level);
			// Compress chunks of 10k records
			int k = 0;
			for(; k + 10000 < blk->n; k += 10000){
				for(int l = k; l < k + 10000; ++l) obuf << blk->rcds[l];
				rec.nc += zcodec.StreamCompress(obuf, obuf2, ostream, twk1_two_t::packed_size * 5000);
				progress->cmps += 10000;
				obuf.reset(); obuf2.reset();
			}

			// Compress residual records
			//std::cerr << "remainder=" << k << "/" << blk->n << std::endl;

			progress->cmps += blk->n - k;
			for(int l = k; l < blk->n; ++l) obuf << blk->rcds[l];
			rec.nc += zcodec.StreamCompress(obuf, obuf2, ostream, twk1_two_t::packed_size * 5000);
			obuf.reset(); obuf2.reset();

			zcodec.StopStreamCompress();
			//std::cerr << "after stop=" << zcodec.outbuf.pos << "/" << zcodec.outbuf.size << std::endl;
			ostream.write((const char*)zcodec.outbuf.dst, zcodec.outbuf.pos);
			rec.nc += zcodec.outbuf.pos;
			rec.n = blk->n * twk1_two_t::packed_size;
			rec.fend = ostream.tellp();
			local_idx.push_back(rec);

			blk->reset();
		}
		ostream.flush();
		std::cerr << utility::timestamp("LOG","THREAD") << "Finished: " << tmp_filename << " with " << f << "-" << t << ". Sorted n=" << tot << " variants with size=" << utility::ToPrettyDiskString((uint64_t)ostream.tellp()) << std::endl;
		ostream.close();
		obuf.clear();
		obuf2.clear();
		delete it; it = nullptr;
		delete blk; blk = nullptr;

		return true;
	}

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
			const uint32_t n_compressed)
	{
		stream = std::ifstream(file, std::ios::binary);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR") << "Could not open \"" << file << "\"..." << std::endl;
			return false;
		}

		stream.seekg(foff);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to seek in file \"" << file << "\"..." << std::endl;
			return false;
		}

		n = 0; it = 0; n_unc_cum = 0;
		b_left = fend - foff;
		it_tot = 0;
		//std::cerr << "data left=" << b_left << ": " << foff << "-" << fend << std::endl;
		n_unc = n_uncompressed; n_cmp = n_compressed;
		n_tot = n_unc / twk1_two_t::packed_size;
		off_start = foff; off_end = fend;
		out.reset(); out2.reset();

		zcodec.InitStreamDecompress();

		return true;
	}

	bool Next(twk1_two_t& rec, const uint32_t b_read = 256000){
		if(it == n){
			if(NextBlock(b_read) == false)
				return false;
		}

		out2 >> rec;
		++it; ++it_tot;
		return true;
	}

	bool NextBlock(const uint32_t b_read = 256000){
		if(b_left == 0)
			return false;

		const uint32_t residual = out2.size() % twk1_two_t::packed_size;

		if(residual){
			//std::cerr << "has residual=" << residual << std::endl;
			std::memmove(out2.data(), out2.data() + (out2.size() - residual), residual);
		}

		//b_read = b_read > b_left ? b_left : b_read;
		const size_t n_read = b_left > b_read ? b_read : b_left;
		n = 0; it = 0;

		n_unc_cum += out2.size();
		if(out.n_chars_ - out.iterator_position_ != 0){
			//std::cerr << "moving=" << out.n_chars_ - out.iterator_position_ << "(" << out.n_chars_ << " and " << out.iterator_position_ << ")" << std::endl;
			std::memmove(out.data(), out.data() + out.iterator_position_, out.n_chars_ - out.iterator_position_);
			out.n_chars_ = out.n_chars_ - out.iterator_position_;
		} else out.n_chars_ = 0;

		out.iterator_position_ = 0;
		out2.n_chars_ = residual;
		out2.iterator_position_ = 0;
		n_unc_cum -= residual;

		if(out.capacity() == 0 || b_read > out.capacity()){
			out.resize(b_read);
			out2.resize(b_read);
		}

		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR") << "Stream is dead!" << std::endl;
			return false;
		}

		// Read a block of data
		//std::cerr << "reader=" << stream.tellg() << " left=" << b_left << std::endl;
		stream.read(out.data() + out.n_chars_, n_read);
		if(stream.good() == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to block of data..." << std::endl;
			return false;
		}
		out.n_chars_ += n_read;

		// Decompress
		//std::cerr << "read=" << n_read << " and " << out.iterator_position_ << "/" << out.n_chars_ << std::endl;
		if(zcodec.StreamDecompress(out, out2) == false){
			std::cerr << utility::timestamp("ERROR") << "Fail stream decompression..." << std::endl;
			return false;
		}

		// If input data is not all consumed
		//std::cerr << "decomp=" << out.iterator_position_ << " -> " << out2.size() << std::endl;

		n = out2.size() / twk1_two_t::packed_size;
		//std::cerr << "n=" << n << " it=" << it << std::endl;

		b_left -= n_read;
		return true;
	}

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

struct two_sorter_settings {
	two_sorter_settings() : memory_limit(0.5), c_level(1), n_threads(std::thread::hardware_concurrency()){}

	std::string in, out;
	float memory_limit;
	int c_level, n_threads;
};

class two_sorter {
public:
	bool Sort(const two_sorter_settings& settings){
		this->settings = settings;
		return(Sort());
	}

	bool Sort(){
		if(settings.in.length() == 0){
			std::cerr << utility::timestamp("ERROR") << "No input value specified..." << std::endl;
			return false;
		}

		// File reader.
		two_reader oreader;
		if(oreader.Open(settings.in) == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to open \"" << settings.in << "\"..."  << std::endl;
			return false;
		}

		twk1_two_block_t blk2;
		twk_buffer_t obuf, obuf2;

		// Distrubution.
		uint64_t b_unc = 0, n_recs = 0;
		std::cerr << utility::timestamp("LOG") << "Blocks: " << utility::ToPrettyString(oreader.index.n) << std::endl;
		for(int i = 0; i < oreader.index.n; ++i){
			b_unc  += oreader.index.ent[i].b_unc;
			n_recs += oreader.index.ent[i].n;
		}
		std::cerr << utility::timestamp("LOG") << "Uncompressed size: " << utility::ToPrettyDiskString(b_unc) << std::endl;
		std::cerr << utility::timestamp("LOG") << "Sorting " << utility::ToPrettyString(n_recs) << " records..." << std::endl;

		if(b_unc == 0){
			std::cerr << utility::timestamp("ERROR") << "Cannot sort empty file..." << std::endl;
			return false;
		}

		if(oreader.index.n < settings.n_threads) settings.n_threads = oreader.index.n;
		uint64_t b_unc_thread = b_unc / settings.n_threads;
		std::cerr << utility::timestamp("LOG","THREAD") << "Data/thread: " << utility::ToPrettyDiskString(b_unc_thread) << std::endl;

		std::vector< std::pair<uint32_t,uint32_t> > ranges;
		uint64_t f = 0, t = 0, b_unc_tot = 0;
		for(int i = 0; i < oreader.index.n; ++i){
			if(b_unc_tot >= b_unc_thread){
				ranges.push_back(std::pair<uint32_t,uint32_t>(f, t));
				b_unc_tot = 0;
				f = t;
			}
			b_unc_tot += oreader.index.ent[i].b_unc;
			++t;
		}
		if(f != t){
			ranges.push_back(std::pair<uint32_t,uint32_t>(f, t));
			b_unc_tot = 0;
			f = t;
		}
		assert(ranges.back().second == oreader.index.n);
		assert(ranges.size() <= settings.n_threads);

		twk_sort_progress progress_sort;
		progress_sort.n_cmps = n_recs;
		std::thread* psthread = progress_sort.Start();

		twk_sort_slave* slaves = new twk_sort_slave[settings.n_threads];
		uint32_t range_thread = oreader.index.n / settings.n_threads;
		for(int i = 0; i < settings.n_threads; ++i){
			slaves[i].f = ranges[i].first;
			slaves[i].t = ranges[i].second;
			slaves[i].m_limit = settings.memory_limit;
			slaves[i].filename = settings.in;
			slaves[i].c_level = settings.c_level;
			slaves[i].progress = &progress_sort;

			std::string suffix    = twk_two_writer_t::RandomSuffix();
			std::string base_path = twk_two_writer_t::GetBasePath(settings.out);
			std::string base_name = twk_two_writer_t::GetBaseName(settings.out);
			std::string temp_out  = (base_path.size() ? base_path + "/" : "") + base_name + "_" + suffix + ".two";
			slaves[i].tmp_filename = temp_out;
			std::cerr << utility::timestamp("LOG","THREAD") << "Slave-" << i << ": range=" << slaves[i].f << "->" << slaves[i].t << "/" << oreader.index.n << " and name " << slaves[i].tmp_filename << std::endl;
		}

		for(int i = 0; i < settings.n_threads; ++i){
			if(slaves[i].Start(oreader) == nullptr){
				std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to spawn slave" << std::endl;
				return false;
			}
		}
		for(int i = 0; i < settings.n_threads; ++i) slaves[i].thread->join();
		progress_sort.is_ticking = false;
		progress_sort.PrintFinal();

		// temp
		for(int i = 0; i < settings.n_threads; ++i){
			std::cerr << i << "\t" << slaves[i].run_ivals.size() << std::endl;
			for(int j = 0; j < slaves[i].run_ivals.size(); ++j){
				std::cerr << "\tblock-" << j << ": " << slaves[i].run_ivals[j].size() << "\t(" << slaves[i].run_ivals[j][0].ref_rid << "," << slaves[i].run_ivals[j][0].n_run << "," << slaves[i].run_ivals[j][0].minp << "-" << slaves[i].run_ivals[j][0].maxp << ")";
				for(int k = 1; k < slaves[i].run_ivals[j].size(); ++k){
					std::cerr << ", (" << slaves[i].run_ivals[j][k].ref_rid << "," << slaves[i].run_ivals[j][k].n_run << "," << slaves[i].run_ivals[j][k].minp << "-" << slaves[i].run_ivals[j][k].maxp << ")";
				}
				std::cerr << std::endl;
			}
		}
		//

		uint32_t n_queues = 0;
		for(int i = 0; i < settings.n_threads; ++i){
			for(int j = 0; j < slaves[i].local_idx.size(); ++j){
				++n_queues;
			}
		}


		// Merge
		obuf.reset(); obuf2.reset();
		obuf.resize(256000);
		obuf2.resize(256000);

		//uint32_t k = 0;
		uint64_t maxmem_queue = settings.memory_limit * settings.n_threads * 1e9 / 15; // assume compression ratio is 15
		uint64_t mem_queue = maxmem_queue / n_queues;
		mem_queue = mem_queue < sizeof(twk1_two_t) ? sizeof(twk1_two_t) : mem_queue;

		std::priority_queue<two_queue_entry> queue;

		std::cerr << utility::timestamp("LOG") << "Spawning " << utility::ToPrettyString(n_queues) << " queues with " << utility::ToPrettyDiskString(mem_queue) << " each..." << std::endl;
		twk_two_stream_iterator* its = new twk_two_stream_iterator[n_queues];
		twk1_two_t rec;
		uint64_t n_rec_total = 0;
		uint32_t local_queue = 0;
		for(int i = 0; i < settings.n_threads; ++i){
			for(int j = 0; j < slaves[i].local_idx.size(); ++j, ++local_queue){
				// open iterators
				if(its[local_queue].Open(slaves[i].tmp_filename,
										 slaves[i].local_idx[j].foff,
										 slaves[i].local_idx[j].fend,
										 slaves[i].local_idx[j].n,
										 slaves[i].local_idx[j].nc) == false)
				{
					std::cerr << utility::timestamp("ERROR") << "Failed open \"" << slaves[i].tmp_filename << "\"..." << std::endl;
					return false;
				}

				if(its[local_queue].Next(rec, mem_queue) == false){
					std::cerr << utility::timestamp("ERROR") << "Failed to get next" << std::endl;
					return false;
				}

				queue.push(two_queue_entry(rec, local_queue));

				n_rec_total += slaves[i].local_idx[j].n / twk1_two_t::packed_size;
			}
		}

		if(queue.empty()){
			std::cerr << utility::timestamp("ERROR","SORT") << "No data in queue..." << std::endl;
			return false;
		}

		twk_two_writer_t owriter;
		owriter.oindex.SetChroms(oreader.hdr.GetNumberContigs());
		if(settings.out.size() == 0 || (settings.out.size() == 1 && settings.out == "-")){
			std::cerr << utility::timestamp("LOG","WRITER") << "Writing to stdout..." << std::endl;
	 	} else {
	 		std::string extension = twk_two_writer_t::GetExtension(settings.out);
			if(extension != "two"){
				settings.out += ".two";
			}
			std::cerr << utility::timestamp("LOG","WRITER") << "Opening \"" << settings.out << "\"..." << std::endl;
	 	}

		if(owriter.Open(settings.out) == false){
			std::cerr << utility::timestamp("ERROR") << "Failed top open \"" << settings.out << "\"..." << std::endl;
			return false;
		}
		owriter.mode = 'b';
		owriter.oindex.state = TWK_IDX_SORTED;
		owriter.SetCompressionLevel(settings.c_level);
		// Write header
		std::string sort_string = "\n##tomahawk_sortVersion=" + std::string(VERSION) + "\n";
		sort_string += "##tomahawk_sortCommand=" + LITERAL_COMMAND_LINE + "; Date=" + utility::datetime() + "\n";
		oreader.hdr.literals_ += sort_string;
		if(owriter.WriteHeader(oreader) == false){
			std::cerr << "failed to write header" << std::endl;
			return false;
		}

		// Reference
		uint32_t ridA = queue.top().rec.ridA;

		Timer timer; timer.Start();
		uint64_t n_entries_out = 0;
		twk_sort_progress progress;
		progress.n_cmps = n_recs;
		std::thread* pthread = progress.Start();

		while(queue.empty() == false){
			// peek at top entry in queue
			const uint32_t id = queue.top().qid;

			if(queue.top().rec.ridA != ridA){
				if(owriter.WriteBlock() == false){
					std::cerr << utility::timestamp("ERROR") << "Failed to flush block..." << std::endl;
					return false;
				}
			}
			owriter.Add(queue.top().rec);
			ridA = queue.top().rec.ridA;

			++progress.cmps;

			// remove this record from the queue
			queue.pop();

			while(its[id].Next(rec, mem_queue)){
				if(!(rec < queue.top().rec)){
					queue.push( two_queue_entry(rec, id) );
					break;
				}

				if(rec.ridA != ridA){
					if(owriter.WriteBlock() == false){
						std::cerr << utility::timestamp("ERROR") << "Failed to flush block..." << std::endl;
						return false;
					}
				}
				owriter.Add(rec);
				ridA = rec.ridA;

				++progress.cmps;
			}
		}
		progress.is_ticking = false;
		progress.PrintFinal();

		owriter.flush();
		owriter.WriteFinal();
		owriter.close();
		std::cerr << utility::timestamp("LOG") << "Finished merging! Time: " << timer.ElapsedString() << std::endl;
		//std::cerr << "deleting intermediary" << std::endl;

		std::cerr << utility::timestamp("LOG") << "Deleting temp files..." << std::endl;
		std::cerr.flush();
		for(int i = 0; i < settings.n_threads; ++i){
		if( remove( slaves[i].tmp_filename.c_str() ) != 0 ){
			std::cerr << utility::timestamp("ERROR") << "Error deleting file " << slaves[i].tmp_filename << "!" << std::endl;
		} else {
			std::cerr << utility::timestamp("LOG") << "Deleted " << slaves[i].tmp_filename << std::endl;
		  }
		}

		delete[] slaves;
		delete[] its;
		std::cerr << utility::timestamp("LOG") << "Finished!" << std::endl;
		return true;
	}

	two_sorter_settings settings;
};

}



#endif /* INCLUDE_TWO_SORTER_H_ */
