/*
Copyright (C) 2016-present Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/
#include <getopt.h>

#include "utility.h"
#include "two_reader.h"

void sort_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Convert binary TWO->LD/TWO, subset and slice TWO data\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " view [options] -i <in.two>\n\n"
	"Options:\n"
	"  -i FILE   input TWO file (required)\n"
	"  -o FILE   output file (- for stdout; default: -)\n"
	"  -m FLOAT  maximum memory usage per thread in GB (default: 0.5)\n\n";
}

int sort(int argc, char** argv){
	if(argc < 3){
		sort_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{"output",      optional_argument, 0, 'o' },
		{"memory-usage", optional_argument, 0, 'm' },


		{0,0,0,0}
	};

	std::string in, out;
	float memory_limit = 0.5;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:o:m:?", long_options, &long_index)) != -1){
		hits += 2;
		switch (c){
		case ':':   /* missing option argument */
			fprintf(stderr, "%s: option `-%c' requires an argument\n",
					argv[0], optopt);
			break;

		case '?':
		default:
			fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
					argv[0], optopt);
			break;

		case 'i':
			in = std::string(optarg);
			break;
		case 'o':
			out = std::string(optarg);
			break;
		case 'm':
			memory_limit = atof(optarg);
			break;


		}
	}

	if(in.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	tomahawk::two_reader oreader;

	if(oreader.Open(in) == false){
		std::cerr << "failed to open" << std::endl;
		return false;
	}

	std::string extension = tomahawk::twk_two_writer_t::GetExtension(out);
	if(extension != "two"){
		out += ".two";
	}
	tomahawk::twk1_two_block_t blk2;
	tomahawk::twk_buffer_t obuf, obuf2;

	struct sort_helper {
		sort_helper() : foff(0), fend(0), n(0), nc(0){}

		uint64_t foff, fend, n, nc;
	};

	struct twk_sort_slave {

		std::thread* Start(tomahawk::two_reader& rdr){
			blk.resize((m_limit*1e9)/sizeof(tomahawk::twk1_two_t));
			if(f > t) return nullptr;
			if(t - f == 0) return nullptr;
			//delete thread;

			stream = std::ifstream(filename,std::ios::binary | std::ios::in);
			if(stream.good() == false){
				std::cerr << "failed to open=" << filename << std::endl;
				return nullptr;
			}

			stream.seekg(rdr.index.ent[f].foff);
			if(stream.good() == false){
				std::cerr << "failed to seek to pos=" << rdr.index.ent[f].foff << " in "  << filename << std::endl;
				return nullptr;
			}

			it.stream = &stream;

			ostream = std::ofstream(tmp_filename, std::ios::binary | std::ios::out);
			if(ostream.good() == false){
				std::cerr << "failed to open temp output file=" << tmp_filename << std::endl;
				return nullptr;
			}

			thread = new std::thread(&twk_sort_slave::Sort, this);
			return(thread);
		}

		bool Sort(){
			uint32_t tot = 0;
			for(int i = f; i < t; ++i){
				assert(it.NextBlock());
				tot += it.GetBlock().n;
				for(int j = 0; j < it.blk.n; ++j){
					//it.blk[j].Print(std::cerr);
					if(blk.n == blk.m){
						blk.Sort();

						sort_helper rec;
						rec.foff = ostream.tellp();
						zcodec.InitStreamCompress(1);
						// Compress chunks of 10k records
						int k = 0;
						for(; k + 10000 < blk.n; k += 10000){
							for(int l = k; l < k + 10000; ++l) obuf << blk.rcds[l];
							rec.nc += zcodec.StreamCompress(obuf, obuf2, ostream, tomahawk::twk1_two_t::packed_size * 5000);
							obuf.reset(); obuf2.reset();
						}

						// Compress residual records
						//std::cerr << "remainder=" << k << "/" << blk.n << std::endl;

						for(int l = k; l < blk.n; ++l) obuf << blk.rcds[l];
						rec.nc += zcodec.StreamCompress(obuf, obuf2, ostream, tomahawk::twk1_two_t::packed_size * 5000);
						obuf.reset(); obuf2.reset();

						zcodec.StopStreamCompress();
						//std::cerr << "after stop=" << zcodec.outbuf.pos << "/" << zcodec.outbuf.size << std::endl;
						ostream.write((const char*)zcodec.outbuf.dst, zcodec.outbuf.pos);
						rec.nc += zcodec.outbuf.pos;
						rec.n = blk.n * tomahawk::twk1_two_t::packed_size;
						rec.fend = ostream.tellp();
						local_idx.push_back(rec);

						blk.reset();

					}
					blk += it.blk[j];
					//it.blk[j].Print(std::cout);
				}
			}
			// any possible remainder
			if(blk.n){
				//std::cerr << "resseting" << std::endl;
				blk.Sort();
				//obuf << blk2; // poor memory management:

				sort_helper rec;
				rec.foff = ostream.tellp();
				zcodec.InitStreamCompress(1);
				// Compress chunks of 10k records
				int k = 0;
				for(; k + 10000 < blk.n; k += 10000){
					for(int l = k; l < k + 10000; ++l) obuf << blk.rcds[l];
					rec.nc += zcodec.StreamCompress(obuf, obuf2, ostream, tomahawk::twk1_two_t::packed_size * 5000);
					obuf.reset(); obuf2.reset();
				}

				// Compress residual records
				//std::cerr << "remainder=" << k << "/" << blk.n << std::endl;

				for(int l = k; l < blk.n; ++l) obuf << blk.rcds[l];
				rec.nc += zcodec.StreamCompress(obuf, obuf2, ostream, tomahawk::twk1_two_t::packed_size * 5000);
				obuf.reset(); obuf2.reset();

				zcodec.StopStreamCompress();
				//std::cerr << "after stop=" << zcodec.outbuf.pos << "/" << zcodec.outbuf.size << std::endl;
				ostream.write((const char*)zcodec.outbuf.dst, zcodec.outbuf.pos);
				rec.nc += zcodec.outbuf.pos;
				rec.n = blk.n * tomahawk::twk1_two_t::packed_size;
				rec.fend = ostream.tellp();
				local_idx.push_back(rec);

				blk.reset();
			}
			ostream.flush();
			std::cerr << "Finished: " << tmp_filename << " with " << f << "-" << t << std::endl;
			std::cerr << "Sorted n=" << tot << " variants with size=" << tomahawk::utility::ToPrettyDiskString((uint64_t)ostream.tellp()) << std::endl;
			ostream.close();

			return true;
		}

		float m_limit;
		uint32_t n; // limit in gb, number of entries that corresponds to
		uint32_t f, t, bl_size; // (from,to)-tuple, flush block-size
		std::string filename, tmp_filename; // temporary filename
		std::thread* thread;
		std::ifstream stream;
		std::ofstream ostream;
		tomahawk::twk1_two_iterator it;
		tomahawk::twk1_two_block_t blk;
		std::vector<sort_helper> local_idx; // local offset index
		tomahawk::ZSTDCodec zcodec;
		tomahawk::twk_buffer_t obuf, obuf2;
	};

	// Distrubution.
	uint64_t b_unc = 0;
	std::cerr << "oreader size=" << oreader.index.n << std::endl;
	for(int i = 0; i < oreader.index.n; ++i){
		b_unc += oreader.index.ent[i].b_unc;
		std::cerr << oreader.index.ent[i].foff << " " << oreader.index.ent[i].b_cmp << " and " << oreader.index.ent[i].b_unc << std::endl;
	}
	std::cerr << "uncompressed size=" << b_unc << std::endl;

	uint32_t n_threads = std::thread::hardware_concurrency();
	if(oreader.index.n < n_threads) n_threads = oreader.index.n;
	uint32_t b_unc_thread = b_unc / n_threads;
	std::cerr << "bytes / thread = " << b_unc_thread << std::endl;
	return(0);
	std::vector< std::pair<uint32_t,uint32_t> > ranges;
	uint32_t f = 0, t = 0, b_unc_tot = 0;
	for(int i = 0; i < oreader.index.n; ++i){
		if(b_unc_tot >= b_unc_thread){
			std::cerr << "adding=" << f << "-" << t <<  " with " << b_unc_tot << std::endl;
			ranges.push_back(std::pair<uint32_t,uint32_t>(f, t));
			b_unc_tot = 0;
			f = t;
		}
		b_unc_tot += oreader.index.ent[i].b_unc;
		++t;
	}
	if(f != t){
		std::cerr << "adding=" << f << "-" << t <<  " with " << b_unc_tot << std::endl;
		ranges.push_back(std::pair<uint32_t,uint32_t>(f, t));
		b_unc_tot = 0;
		f = t;
	}
	std::cerr << "ranges=" << ranges.size() << std::endl;

	twk_sort_slave* slaves = new twk_sort_slave[n_threads];
	std::cerr << "index=" << oreader.index.n << " -> " << oreader.index.n / n_threads << std::endl;
	uint32_t range_thread = oreader.index.n / n_threads;
	for(int i = 0; i < n_threads; ++i){
		slaves[i].f = range_thread * i;
		slaves[i].t = (i + 1 == n_threads ? oreader.index.n : (range_thread * (i+1)));
		slaves[i].m_limit = memory_limit;
		slaves[i].filename = in;

		std::string suffix    = tomahawk::twk_two_writer_t::RandomSuffix();
		std::string base_path = tomahawk::twk_two_writer_t::GetBasePath(out);
		std::string base_name = tomahawk::twk_two_writer_t::GetBaseName(out);
		std::string temp_out  = (base_path.size() ? base_path + "/" : "") + base_name + "_" + suffix + ".two";
		slaves[i].tmp_filename = temp_out;
		std::cerr << "Slave-" << i << ": range=" << slaves[i].f << "->" << slaves[i].t << "/" << oreader.index.n << " and name " << slaves[i].tmp_filename << std::endl;
	}

	for(int i = 0; i < n_threads; ++i){
		if(slaves[i].Start(oreader) == nullptr){
			std::cerr << "failed to spawn slave" << std::endl;
			return 1;
		}
	}
	for(int i = 0; i < n_threads; ++i) slaves[i].thread->join();
	for(int i = 0; i < n_threads; ++i){
		slaves[i].obuf.clear(); // release memory
		slaves[i].obuf2.clear();
		slaves[i].it.buf.clear();
		slaves[i].it.oblk.bytes.clear();
	}


	uint32_t n_queues = 0;
	for(int i = 0; i < n_threads; ++i){
		for(int j = 0; j < slaves[i].local_idx.size(); ++j){
			std::cerr << i << ": " << slaves[i].local_idx[j].foff << "->" << slaves[i].local_idx[j].fend << " " << slaves[i].local_idx[j].n << " " << slaves[i].local_idx[j].nc << std::endl;
			++n_queues;
		}
	}

	struct twk_two_stream_iterator {

		bool Open(const std::string file,
				const uint64_t foff,
				const uint64_t fend,
				const uint32_t n_uncompressed,
				const uint32_t n_compressed)
		{
			stream = std::ifstream(file, std::ios::binary);
			if(stream.good() == false){
				std::cerr << "could not open " << file << std::endl;
				return false;
			}

			stream.seekg(foff);
			if(stream.good() == false){
				std::cerr << "seek failed" << std::endl;
			}

			n = 0; it = 0; n_unc_cum = 0;
			b_left = fend - foff;
			it_tot = 0;
			//std::cerr << "data left=" << b_left << ": " << foff << "-" << fend << std::endl;
			n_unc = n_uncompressed; n_cmp = n_compressed;
			n_tot = n_unc / tomahawk::twk1_two_t::packed_size;
			off_start = foff; off_end = fend;
			out.reset(); out2.reset();

			zcodec.InitStreamDecompress();

			return true;
		}

		bool Next(tomahawk::twk1_two_t& rec){
			if(it == n){
				if(NextBlock() == false)
					return false;
			}

			out2 >> rec;
			++it; ++it_tot;
			return true;
		}

		bool NextBlock(uint32_t b_read = 256000){
			if(b_left == 0)
				return false;

			const uint32_t residual = out2.size() % tomahawk::twk1_two_t::packed_size;

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
				std::cerr << "stream is dead" << std::endl;
				return false;
			}

			// Read a block of data
			//std::cerr << "reader=" << stream.tellg() << " left=" << b_left << std::endl;
			stream.read(out.data() + out.n_chars_, n_read);
			if(stream.good() == false){
				std::cerr << "failed to read: " << n_read << std::endl;
				return false;
			}
			out.n_chars_ += n_read;

			// Decompress
			//std::cerr << "read=" << n_read << " and " << out.iterator_position_ << "/" << out.n_chars_ << std::endl;
			if(zcodec.StreamDecompress(out, out2) == false){
				std::cerr << "fail stream decompress" << std::endl;
				return false;
			}

			// If input data is not all consumed
			//std::cerr << "decomp=" << out.iterator_position_ << " -> " << out2.size() << std::endl;

			n = out2.size() / tomahawk::twk1_two_t::packed_size;
			//std::cerr << "n=" << n << " it=" << it << std::endl;

			b_left -= n_read;
			return true;
		}

		uint32_t n, it, n_tot, it_tot; // n_records, it_pos
		uint32_t b_left, n_unc, n_cmp, n_unc_cum;
		uint64_t off_start, off_end;
		std::ifstream stream;
		tomahawk::twk_buffer_t out, out2;
		tomahawk::ZSTDCodec zcodec;
	};

	struct two_queue_entry {
	public:
		two_queue_entry(const tomahawk::twk1_two_t& data, uint32_t streamID) :
			qid(streamID), rec(data)
	    {}

	    inline bool operator<(const two_queue_entry& other) const {
	        return(!(rec < other.rec));
	    }

	public:
	    uint32_t qid;
	    tomahawk::twk1_two_t rec;
	};

	obuf.reset(); obuf2.reset();
	obuf.resize(256000);
	obuf2.resize(256000);

	//uint32_t k = 0;
	std::priority_queue<two_queue_entry> queue;

	std::cerr << "spawning=" << n_queues << " queues" << std::endl;
	twk_two_stream_iterator* its = new twk_two_stream_iterator[n_queues];
	tomahawk::twk1_two_t rec;
	uint64_t n_rec_total = 0;
	uint32_t local_queue = 0;
	for(int i = 0; i < n_threads; ++i){
		for(int j = 0; j < slaves[i].local_idx.size(); ++j, ++local_queue){
			// open iterators
			if(its[local_queue].Open(slaves[i].tmp_filename, slaves[i].local_idx[j].foff, slaves[i].local_idx[j].fend, slaves[i].local_idx[j].n, slaves[i].local_idx[j].nc) == false){
				std::cerr << "failed open" << std::endl;
				return false;
			}

			if(its[local_queue].Next(rec) == false){
				std::cerr << "failed to get next" << std::endl;
				return false;
			}
			std::cerr << "queue-" << local_queue << std::endl;
			//rec.PrintLD(std::cerr);
			queue.push(two_queue_entry(rec, local_queue));

			n_rec_total += slaves[i].local_idx[j].n / tomahawk::twk1_two_t::packed_size;
			std::cerr << "total=" << its[local_queue].n_tot << std::endl;
		}
	}

	if(queue.empty()){
		std::cerr << tomahawk::utility::timestamp("ERROR","SORT") << "No data in queue..." << std::endl;
		return false;
	}

	//delete[] slaves;
	//return 0;

	//std::cerr << "top" << std::endl;
	//queue.top().rec.PrintLD(std::cerr);

	tomahawk::twk_two_writer_t owriter;
	owriter.oindex.SetChroms(oreader.hdr.GetNumberContigs());
	std::cerr << "opening=" << out << std::endl;
	if(owriter.Open(out) == false){
		std::cerr << "failed top open" << std::endl;
		return false;
	}
	owriter.mode = 'b';
	owriter.oindex.state = TWK_IDX_SORTED;
	// Write header
	std::string sort_string = "\n##tomahawk_sortVersion=" + std::to_string(VERSION) + "\n";
	sort_string += "##tomahawk_sortCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + tomahawk::utility::datetime() + "\n";
	oreader.hdr.literals_ += sort_string;
	if(owriter.WriteHeader(oreader) == false){
		std::cerr << "failed to write header" << std::endl;
		return false;
	}


	// Reference
	uint32_t ridA = queue.top().rec.ridA;

	tomahawk::Timer timer; timer.Start();
	uint64_t n_entries_out = 0;
	while(queue.empty() == false){
		// peek at top entry in queue
		const uint32_t id = queue.top().qid;

		if(queue.top().rec.ridA != ridA){
			//std::cerr << "change in ridA=" << queue.top().rec.ridA << "!=" << ridA << " -> " << owriter.oblock.n << std::endl;
			if(owriter.WriteBlock() == false){
				std::cerr << "failed to flush" << std::endl;
				return false;
			}
		}

		//queue.top().rec.PrintLD(std::cout);
		owriter.Add(queue.top().rec);
		ridA = queue.top().rec.ridA;


		++n_entries_out;
		if(n_entries_out % 1000000 == 0)
			std::cerr << tomahawk::utility::timestamp("LOG") << tomahawk::utility::ToPrettyString(n_entries_out) << "/" << tomahawk::utility::ToPrettyString(n_rec_total) << std::endl;

		// remove this record from the queue
		queue.pop();

		while(its[id].Next(rec)){
			if(!(rec < queue.top().rec)){
				queue.push( two_queue_entry(rec, id) );
				break;
			}

			if(rec.ridA != ridA){
				//std::cerr << "change in ridA=" << rec.ridA << "!=" << ridA << " -> " << owriter.oblock.n << std::endl;

				if(owriter.WriteBlock() == false){
					std::cerr << "failed to flush" << std::endl;
					return false;
				}
			}
			//rec.PrintLD(std::cout);
			owriter.Add(rec);
			ridA = rec.ridA;

			++n_entries_out;
			if(n_entries_out % 1000000 == 0)
				std::cerr << tomahawk::utility::timestamp("LOG") << tomahawk::utility::ToPrettyString(n_entries_out) << "/" << tomahawk::utility::ToPrettyString(n_rec_total) << std::endl;
		}
	}

	owriter.flush();
	owriter.WriteFinal();
	owriter.close();
	std::cerr << "merge time=" << timer.ElapsedString() << std::endl;
	std::cerr << "deleting intermediary" << std::endl;

	for(int i = 0; i < n_threads; ++i){
	if( remove( slaves[i].tmp_filename.c_str() ) != 0 ){
	    std::cerr << "Error deleting file=" << slaves[i].tmp_filename << std::endl;
	    //return 1;
	} else {
	    std::cerr << "File successfully deleted" << slaves[i].tmp_filename << std::endl;
	  }
	}

	delete[] slaves;
	delete[] its;
	return 0;
}
