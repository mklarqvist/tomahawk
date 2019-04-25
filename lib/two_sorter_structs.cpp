#include "two_sorter_structs.h"
#include "two_reader.h"

namespace tomahawk {

std::thread* twk_sort_slave::Start(IndexOutput& rdr){
	if(f > t) return nullptr;
	if(t - f == 0) return nullptr;
	//delete thread;

	stream.open(filename,std::ios::binary | std::ios::in);
	if(stream.good() == false){
		std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to open \"" << filename << "\"..." << std::endl;
		return nullptr;
	}

	stream.seekg(rdr.ent[f].foff);
	if(stream.good() == false){
		std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to seek to position " << rdr.ent[f].foff << " in \""  << filename << "\"..." << std::endl;
		return nullptr;
	}

	it  = new twk1_two_iterator;
	blk = new twk1_two_block_t;
	blk->resize((m_limit*1e9)/sizeof(twk1_two_t));
	it->stream = &stream;

	ostream.open(tmp_filename, std::ios::binary | std::ios::out);
	if(ostream.good() == false){
		std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to open temp output file \"" << tmp_filename << "\"..." << std::endl;
		return nullptr;
	}

	thread = new std::thread(&twk_sort_slave::Sort, this);
	return(thread);
}

bool twk_sort_slave::Sort(){
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
				zcodec.WriteOutbuf(ostream);
				//ostream.write((const char*)zcodec.outbuf.dst, zcodec.outbuf.pos);
				rec.nc += zcodec.GetOutputSize();
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
		zcodec.WriteOutbuf(ostream);
		rec.nc += zcodec.GetOutputSize();
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



bool twk_two_stream_iterator::Open(const std::string file,
		const uint64_t foff,
		const uint64_t fend,
		const uint32_t n_uncompressed,
		const uint32_t n_compressed)
{
	stream.open(file, std::ios::binary);
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

bool twk_two_stream_iterator::Next(twk1_two_t& rec, const uint32_t b_read){
	if(it == n){
		if(NextBlock(b_read) == false)
			return false;
	}

	out2 >> rec;
	++it; ++it_tot;
	return true;
}

bool twk_two_stream_iterator::NextBlock(const uint32_t b_read){
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

}
