#ifndef IMPORTER_H_
#define IMPORTER_H_

#include <omp.h>

#include "zstd.h"

#include "tomahawk.h"
#include "vcf_reader.h"
#include "generic_iterator.h"
#include "buffer.h"
#include "genotype_encoder.h"
#include "core.h"
#include "index.h"
#include "zstd_codec.h"

// temp
#include "twk_reader.h"
#include "ld.h"
#include "timer.h"
#include "writer.h"
#include "core.h"

namespace tomahawk {

struct twk_vimport_settings {
	twk_vimport_settings() : block_size(500), input("-"), output("-"){}

	uint32_t block_size;
	std::string input, output;
};

class twk_variant_importer {
public:
	bool Import(twk_vimport_settings& settings){
		this->settings = settings;
		return(this->Import());
	}

	bool Import(void){
		// Start timer.
		Timer timer; timer.Start();

		// Retrieve a unique VcfReader.
		std::string filename = "/home/mk21/Downloads/1kgp3_chr20.bcf";
		//std::string filename = "/media/mdrk/NVMe/1kgp3/bcf/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf";
		//filename = "/home/mk21/Downloads/fish_callset_filtered.bcf";
		//filename = "/home/mk21/Downloads/randomized.100k-50k.bcf";

		ProgramMessage();
		std::cerr << utility::timestamp("LOG") << "Calling import..." << std::endl;
		std::cerr << utility::timestamp("LOG") << "Opening " << filename << "..." << std::endl;
		std::unique_ptr<tomahawk::io::VcfReader> vcf = tomahawk::io::VcfReader::FromFile(filename, 8);
		if(vcf == nullptr){
			std::cerr << "failed to get vcfreader" << std::endl;
			return false;
		}

		if(vcf->vcf_header_.GetFormat("GT") == nullptr){
			std::cerr << "Genotype data not set in this file" << std::endl;
			return false;
		}
		std::cerr << utility::timestamp("LOG","VCF") << "Constructing lookup table for " << utility::ToPrettyString(vcf->vcf_header_.GetNumberContigs()) << " contigs..." << std::endl;
		std::cerr << utility::timestamp("LOG","VCF") << "Samples: " << utility::ToPrettyString(vcf->vcf_header_.GetNumberSamples()) << "..." << std::endl;

		// Todo: add header literal tracing input parameters
		tomahawk::Index index(vcf->vcf_header_.GetNumberContigs());

		std::ostream* stream = nullptr; bool stream_delete = true;
		std::string outfile = "/home/mk21/Downloads/debug.twk";
		outfile = "-";
		if(outfile.size() == 0 || (outfile.size() == 1 && outfile[0] == '-')){
			stream = &std::cout;
			stream_delete = false;
		}
		else {
			std::cerr << utility::timestamp("LOG","WRITER") << "Opening " << outfile << "..." << std::endl;
			stream = new std::ofstream;
			std::ofstream* outstream = reinterpret_cast<std::ofstream*>(stream);
			outstream->open(outfile,std::ios::out | std::ios::binary);
			if(!outstream->good()){
				std::cerr << "failed to open" << std::endl;
				return false;
			}
		}

		tomahawk::ZSTDCodec zcodec;
		stream->write(tomahawk::TOMAHAWK_MAGIC_HEADER.data(), tomahawk::TOMAHAWK_MAGIC_HEADER_LENGTH);
		tomahawk::twk_buffer_t buf(256000), obuf(256000);
		buf << vcf->vcf_header_;
		//std::cerr << "header buf size =" << buf.size() << std::endl;
		if(zcodec.Compress(buf, obuf, 10) == false){
			std::cerr << "failed to compress" << std::endl;
			return false;
		}
		//std::cerr << buf.size() << "->" << obuf.size() << " -> " << (float)buf.size()/obuf.size() << std::endl;

		stream->write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		stream->write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		stream->write(obuf.data(),obuf.size());
		buf.reset();

		// start
		stream->flush();
		//const uint64_t start_pos = outstream->tellp();

		tomahawk::twk1_block_t block;
		uint32_t n_vnt_dropped = 0;
		uint32_t n_total_rec = 0;

		while(vcf->next(BCF_UN_ALL)){
			tomahawk::twk1_t entry;
			entry.pos = vcf->bcf1_->pos;
			//if(n_total_rec == 200000){
				//std::cerr << "break at 10k" << std::endl;
				//break;
			//}

			if(vcf->bcf1_->n_fmt){
				if(vcf->vcf_header_.GetFormat(vcf->bcf1_->d.fmt[0].id)->id == "GT"){
					if(vcf->bcf1_->d.fmt[0].n != 2){
						std::cerr << "do not support non-diploid" << std::endl;
						++n_vnt_dropped;
						continue;
					}

					if(vcf->bcf1_->n_allele != 2){
						//std::cerr << "allele != 2" << std::endl;
						++n_vnt_dropped;
						continue;
					}

					bool valid = true;
					for(int i = 0; i < 2; ++i){
						if(std::regex_match(vcf->bcf1_->d.allele[i],tomahawk::TWK_REGEX_CANONICAL_BASES) == false){
							//std::cerr << "failed match= " << vcf->bcf1_->d.allele[i] << std::endl;
							++n_vnt_dropped;
							valid = false;
							break;
						}
					}
					if(valid == false) continue;

					// Encode alleles.
					entry.EncodeAlleles(vcf->bcf1_->d.allele[0][0],vcf->bcf1_->d.allele[1][0]);
					assert(entry.GetAlleleA() == vcf->bcf1_->d.allele[0][0]);
					assert(entry.GetAlleleB() == vcf->bcf1_->d.allele[1][0]);

					// Encode genotypes.
					if(tomahawk::GenotypeEncoder::Encode(vcf->bcf1_, entry) == false){
						std::cerr << "invalid encoding" << std::endl;
						++n_vnt_dropped;
						continue;
					}
					entry.rid = vcf->bcf1_->rid;

					if(block.n != 0){
						if(block.rid != vcf->bcf1_->rid){
							std::cerr << "wrong rid@" << block.rid << "!=" << vcf->bcf1_->rid << std::endl;
							buf << block;
							/*
							std::cerr << "buf now=" << buf.size() << std::endl;
							twk1_block_t block2;
							buf >> block2;
							std::cerr << "block2=" << block2.n << "==" << block.n << std::endl;
							assert(block2.n == block.n);
							assert(block.size() == block2.size());
							*/
							tomahawk::IndexEntry ent;
							ent.n = block.n; ent.minpos = block.minpos; ent.maxpos = block.maxpos;
							ent.rid = block.rid;
							ent.foff = stream->tellp();
							block.clear();
							block.rid = vcf->bcf1_->rid;
							block.minpos = vcf->bcf1_->pos;

							tomahawk::twk_buffer_t obuf(buf.size() + 65536);
							assert(zcodec.Compress(buf, obuf, 10));
							std::cerr << buf.size() << "->" << obuf.size() << " = " << (float)buf.size()/obuf.size() << std::endl;
							// Write block: size-un, size, buffer
							tomahawk::twk_oblock_t oblock;
							oblock.Write(*stream, buf.size(), obuf.size(), obuf);

							ent.fend = stream->tellp();
							index += ent;
							//std::cerr << tomahawk::utility::timestamp("DEBUG") << outstream->tellp() << std::endl;
							buf.reset();
						}

						if(block.n == 500){
							//std::cerr<< "full=" << block.n << "=" << 500 << std::endl;
							buf << block;
							/*
							std::cerr << "buf now=" << buf.size() << std::endl;
							twk1_block_t block2;
							buf >> block2;
							std::cerr << "block2=" << block2.n << "==" << block.n << std::endl;
							assert(block2.n == block.n);
							assert(block.size() == block2.size());
							*/
							tomahawk::IndexEntry ent;
							ent.n = block.n; ent.minpos = block.minpos; ent.maxpos = block.maxpos;
							ent.rid = block.rid;
							ent.foff = stream->tellp();
							block.clear();
							block.rid = vcf->bcf1_->rid;
							block.minpos = vcf->bcf1_->pos;

							tomahawk::twk_buffer_t obuf(buf.size() + 65536);
							assert(zcodec.Compress(buf, obuf, 10));
							//std::cerr << buf.size() << "->" << obuf.size() << " = " << (float)buf.size()/obuf.size() << std::endl;
							// Write block: size-un, size, buffer
							tomahawk::twk_oblock_t oblock;
							oblock.Write(*stream, buf.size(), obuf.size(), obuf);

							ent.fend = stream->tellp();
							index += ent;

							//std::cerr << tomahawk::utility::timestamp("DEBUG") << outstream->tellp() << std::endl;
							buf.reset();
						}
						++n_total_rec;
						block += entry;
					} else {
						block += entry;
						block.rid = vcf->bcf1_->rid;
						block.minpos = vcf->bcf1_->pos;
						++n_total_rec;
					}
					//entry.gt->print();
				} else {
					++n_vnt_dropped;
					std::cerr << "no genotypes" << std::endl;
				}
			} else {
				++n_vnt_dropped;
				std::cerr << "no fmt" << std::endl;
			}
		}

		// Add last
		if(block.n){
			buf << block;
			tomahawk::IndexEntry ent;
			ent.n = block.n; ent.minpos = block.minpos; ent.maxpos = block.maxpos;
			ent.rid = block.rid;
			ent.foff = stream->tellp();
			block.clear();
			block.rid = vcf->bcf1_->rid;
			block.minpos = vcf->bcf1_->pos;

			tomahawk::twk_buffer_t obuf(buf.size() + 65536);
			assert(zcodec.Compress(buf, obuf, 10));
			//std::cerr << buf.size() << "->" << obuf.size() << " = " << (float)buf.size()/obuf.size() << std::endl;
			// Write block: size-un, size, buffer
			tomahawk::twk_oblock_t oblock;
			oblock.Write(*stream, buf.size(), obuf.size(), obuf);

			ent.fend = stream->tellp();
			index += ent;

			//std::cerr << tomahawk::utility::timestamp("DEBUG") << "Final=" << outstream->tellp() << std::endl;
			buf.reset();
		}

		buf << index;
		//std::cerr << "index buf size =" << buf.size() << std::endl;
		if(zcodec.Compress(buf, obuf, 10) == false){
			std::cerr << "failed to compress" << std::endl;
			return false;
		}
		//std::cerr << buf.size() << "->" << obuf.size() << " -> " << (float)buf.size()/obuf.size() << std::endl;
		const uint64_t offset_start_index = stream->tellp();
		uint8_t marker = 0;
		stream->write(reinterpret_cast<const char*>(&marker),sizeof(uint8_t));
		stream->write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		stream->write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		stream->write(obuf.data(),obuf.size());
		stream->write(reinterpret_cast<const char*>(&offset_start_index),sizeof(uint64_t));
		stream->write(tomahawk::TOMAHAWK_FILE_EOF.data(), tomahawk::TOMAHAWK_FILE_EOF_LENGTH);
		stream->flush();

		std::cerr << utility::timestamp("LOG") << "Wrote: " << utility::ToPrettyString(index.GetTotalVariants()) << " variants to " << utility::ToPrettyString(index.n) << " blocks..." << std::endl;
		std::cerr << utility::timestamp("LOG") << "Finished: " << timer.ElapsedString() << std::endl;

		if(stream_delete) delete stream;
		return(this->Compute());
	}

	bool Compute(){
		std::string filename = "/home/mk21/Downloads/debug.twk";
		std::string outname = "/home/mk21/Downloads/debug.two";

		ProgramMessage();
		std::cerr << utility::timestamp("LOG") << "Calling calc..." << std::endl;
		std::cerr << utility::timestamp("LOG") << "Opening " << filename << "..." << std::endl;

		//*//////////////// Reopen and compute
		twk_reader reader;
		if(reader.Open(filename) == false){
			std::cerr << "failed" << std::endl;
			return false;
		}

		std::cerr << utility::timestamp("LOG") << "Opening " << outname << "..." << std::endl;

		twk_writer_t* writer = new twk_writer_file;
		if(writer->Open(outname) == false){
			std::cerr << "failed to open" << std::endl;
			delete writer;
			return false;
		}

		std::cerr << utility::timestamp("LOG") << "Samples: " << utility::ToPrettyString(reader.hdr.GetNumberSamples()) << std::endl;

		tomahawk::ZSTDCodec zcodec;
		writer->stream->write(tomahawk::TOMAHAWK_LD_MAGIC_HEADER.data(), tomahawk::TOMAHAWK_LD_MAGIC_HEADER_LENGTH);
		tomahawk::twk_buffer_t buf(256000), obuf(256000);
		buf << reader.hdr;
		//std::cerr << "header buf size =" << buf.size() << std::endl;
		assert(zcodec.Compress(buf, obuf, 10));
		//std::cerr << buf.size() << "->" << obuf.size() << " -> " << (float)buf.size()/obuf.size() << std::endl;

		writer->stream->write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		writer->stream->write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		writer->stream->write(obuf.data(),obuf.size());
		writer->stream->flush();
		buf.reset();
		//
		// New index
		IndexOutput index(reader.hdr.GetNumberContigs());

		twk1_blk_iterator bit;
		bit.stream = reader.rstream;

		twk_ld_balancer balancer;
		balancer.Build(reader.index.n, 1, 0);

		Timer timer; timer.Start();
		uint32_t n_blocks = balancer.to - balancer.from;
		twk1_ldd_blk* ldd = new twk1_ldd_blk[n_blocks];
		uint32_t n_variants = 0;
		for(int i = balancer.from; i < balancer.to; ++i) n_variants += reader.index.ent[i].n;
		std::cerr << utility::timestamp("LOG","PARSE") << utility::ToPrettyString(n_variants) << " variants from " << n_blocks << " blocks..." << std::endl;

		bool pre_build = true;
		timer.Start();

#if SIMD_AVAILABLE == 1
		std::cerr << utility::timestamp("LOG","SIMD") << "Vectorized instructions available: " << TWK_SIMD_MAPPING[SIMD_VERSION] << "..." << std::endl;
#else
		std::cerr << utility::timestamp("LOG","SIMD") << "No vectorized instructions available..." << std::endl;
#endif
		std::cerr << utility::timestamp("LOG") << "Constructing list, vector, RLE... ";

		//uint32_t i = 0;
		if(pre_build){
			bit.stream->seekg(reader.index.ent[balancer.from].foff);
			for(int i = 0; i < n_blocks; ++i){
				if(bit.NextBlockRaw() == false){
					std::cerr << "failed to get->" << i << std::endl;
					return false;
				}

				//std::cerr << "next-block=" << bit.stream->tellg() << " n=" << bit.blk.n << std::endl;
				ldd[i].SetOwn(bit, reader.hdr.GetNumberSamples());
				ldd[i].Inflate(reader.hdr.GetNumberSamples(),TWK_LDD_ALL,true);
			}
		} else {
			bit.stream->seekg(reader.index.ent[balancer.from].foff);
			for(int i = 0; i < n_blocks; ++i){
				if(bit.NextBlockRaw() == false){
					std::cerr << "failed to get->" << i << std::endl;
					return false;
				}

				//std::cerr << "next-block=" << bit.stream->tellg() << " n=" << bit.blk.n << std::endl;
				ldd[i].SetOwn(bit, reader.hdr.GetNumberSamples());
			}
		}
		std::cerr << "Done! " << timer.ElapsedString() << std::endl;
		std::cerr << utility::timestamp("LOG","NOTICE") << "Running in fast mode! No matrices will be built..." << std::endl;
		std::cerr << "balancing=" << balancer.from << "-" << balancer.to << std::endl;
		uint64_t n_vnt_cmps = ((uint64_t)n_variants * n_variants - n_variants) / 2;
		std::cerr << utility::timestamp("LOG") << "Performing: " << utility::ToPrettyString(n_vnt_cmps) << " variant comparisons..." << std::endl;

		twk_ld_ticker ticker;
		ticker.i = balancer.from;
		ticker.j = balancer.from;
		ticker.i_start = balancer.from;
		ticker.j_start = balancer.from;
		ticker.limA = balancer.to;
		ticker.limB = balancer.to;
		ticker.ldd  = ldd;
		twk_ld_progress progress;
		progress.n_s = reader.hdr.GetNumberSamples();
		uint32_t n_threads = std::thread::hardware_concurrency();
		//n_threads  = 1;
		twk_ld_slave* slaves = new twk_ld_slave[n_threads];
		std::vector<std::thread*> threads(n_threads);

		std::cerr << utility::timestamp("LOG","THREAD") << "Spawning " << n_threads << " threads: ";
		for(int i = 0; i < n_threads; ++i){
			slaves[i].ldd = ldd;
			slaves[i].n_s = reader.hdr.GetNumberSamples();
			slaves[i].ticker = &ticker;
			slaves[i].engine.SetSamples(reader.hdr.GetNumberSamples());
			slaves[i].engine.SetBlocksize(10000);
			slaves[i].engine.progress = &progress;
			slaves[i].engine.writer   = writer;
			slaves[i].engine.index    = &index;
			slaves[i].progress = &progress;
			threads[i] = slaves[i].Start();
			std::cerr << ".";
		}
		std::cerr << std::endl;

		progress.Start();

		for(int i = 0; i < n_threads; ++i) threads[i]->join();
		for(int i = 0; i < n_threads; ++i) slaves[i].engine.CompressBlock();
		progress.is_ticking = false;
		progress.PrintFinal();
		writer->stream->flush();

		std::cerr << "performed=" << ticker.n_perf << std::endl;


		buf.reset(); obuf.reset();
		buf << index;
		//std::cerr << "index buf size =" << buf.size() << std::endl;
		if(zcodec.Compress(buf, obuf, 10) == false){
			std::cerr << "failed compression" << std::endl;
			return false;
		}
		//std::cerr << buf.size() << "->" << obuf.size() << " -> " << (float)buf.size()/obuf.size() << std::endl;

		const uint64_t offset_start_index = writer->stream->tellp();
		uint8_t marker = 0;
		writer->stream->write(reinterpret_cast<const char*>(&marker),sizeof(uint8_t));
		writer->stream->write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		writer->stream->write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		writer->stream->write(obuf.data(),obuf.size());
		writer->stream->write(reinterpret_cast<const char*>(&offset_start_index),sizeof(uint64_t));
		writer->stream->write(tomahawk::TOMAHAWK_FILE_EOF.data(), tomahawk::TOMAHAWK_FILE_EOF_LENGTH);
		writer->stream->flush();

		delete[] ldd;
		delete[] slaves;
		delete writer;

		return true;

		twk1_ldd_blk blocks[2];

		LDEngine::ep f;

		f[0] = &LDEngine::PhasedVectorized;
		f[1] = &LDEngine::PhasedVectorizedNoMissing;
		f[2] = &LDEngine::PhasedVectorizedNoMissingNoTable;
		f[3] = &LDEngine::UnphasedVectorized;
		f[4] = &LDEngine::UnphasedVectorizedNoMissing;
		f[5] = &LDEngine::PhasedRunlength;
		f[6] = &LDEngine::PhasedList;
		f[7] = &LDEngine::UnphasedRunlength;
		f[8] = &LDEngine::HybridPhased;
		f[9] = &LDEngine::HybridUnphased;

		f[0] = &LDEngine::PhasedList;

		uint8_t unpack_list[10];
		memset(unpack_list, TWK_LDD_VEC, 10);
		unpack_list[5] = TWK_LDD_NONE;
		unpack_list[6] = TWK_LDD_LIST;
		unpack_list[7] = TWK_LDD_NONE;
		unpack_list[8] = TWK_LDD_ALL;
		unpack_list[9] = TWK_LDD_ALL;

		unpack_list[0] = TWK_LDD_LIST;

		ld_perf perfs[10];
		uint32_t perf_size = reader.hdr.GetNumberSamples()*4 + 2;
		for(int i = 0; i < 10; ++i){
			perfs[i].cycles = new uint64_t[perf_size];
			perfs[i].freq   = new uint64_t[perf_size];
			memset(perfs[i].cycles, 0, perf_size*sizeof(uint64_t));
			memset(perfs[i].freq, 0, perf_size*sizeof(uint64_t));
		}

		LDEngine engine;
		engine.SetSamples(reader.hdr.GetNumberSamples());
		engine.SetBlocksize(10000);

		uint32_t dist = 1000000;

		uint32_t from, to; uint8_t type;
		uint64_t n_total = 0;
		uint32_t fl_lim = 0;

		while(true){
			if(!ticker.GetBlockPair(from, to, type)) break;
			if(from == to) std::cerr << utility::timestamp("DEBUG") << from << "," << to << "," << (int)type << std::endl;

			blocks[0] = ldd[from];
			blocks[0].Inflate(reader.hdr.GetNumberSamples(), unpack_list[0]);
			blocks[1] = ldd[to];
			blocks[1].Inflate(reader.hdr.GetNumberSamples(), unpack_list[0]);

			if(type == 1){
				for(int i = 0; i < ldd[from].n_rec; ++i){
					for(int j = i+1; j < ldd[to].n_rec; ++j){
						if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac < 5){
							continue;
						}

						/*
						engine.PhasedVectorized(blocks[0], i, blocks[1], j);
						engine.PhasedVectorizedNoMissing(blocks[0], i, blocks[1], j);
						engine.PhasedVectorizedNoMissingNoTable(blocks[0], i, blocks[1], j);
						engine.UnphasedVectorized(blocks[0], i, blocks[1], j);
						engine.UnphasedVectorizedNoMissing(blocks[0], i, blocks[1], j);
						engine.Runlength(blocks[0], i, blocks[1], j);
						*/
						//engine.AllAlgorithms(blocks[0],i,blocks[0],j,&perfs[method]);
						//(engine.*f[0])(blocks[0], i, blocks[0], j, &perfs[0]);
						engine.PhasedList(blocks[0],i,blocks[0],j,nullptr);
					}
				}
				n_total += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2
			} else {
				for(int i = 0; i < blocks[0].n_rec; ++i){
					for(int j = 0; j < blocks[1].n_rec; ++j){
						if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac < 5 ){
							continue;
						}

						/*
						engine.PhasedVectorized(blocks[0], i, blocks[1], j);
						engine.PhasedVectorizedNoMissing(blocks[0], i, blocks[1], j);
						engine.PhasedVectorizedNoMissingNoTable(blocks[0], i, blocks[1], j);
						engine.UnphasedVectorized(blocks[0], i, blocks[1], j);
						engine.UnphasedVectorizedNoMissing(blocks[0], i, blocks[1], j);
						engine.Runlength(blocks[0], i, blocks[1], j);
						*/
						//(engine.*f[method])(blocks[0], i, blocks[1], j, &perfs[method]);
						engine.PhasedList(blocks[0],i,blocks[1],j,nullptr);
						//engine.AllAlgorithms(blocks[0],i,blocks[1],j,&perfs[method]);
					}
				}
				n_total += blocks[0].n_rec * blocks[1].n_rec;


				if(fl_lim++ == 10){
					std::cerr << "m=" << 666 << " " << utility::ToPrettyString(n_total) << " " << utility::ToPrettyString((uint64_t)((float)n_total/timer.Elapsed().count())) << "/s " << utility::ToPrettyString((uint64_t)((float)n_total*reader.hdr.GetNumberSamples()/timer.Elapsed().count())) << "/gt/s " << timer.ElapsedString() << " " << utility::ToPrettyString(engine.t_out) << std::endl;
					fl_lim = 0;
				}
			}

		}

		std::cerr << "all done" << std::endl;

		//for(int method = 0; method < 10; ++method){
		int method = 0;
			timer.Start();
			n_total = 0; fl_lim = 0;
			for(int b1 = 0; b1 < reader.index.n; ++b1){
				//ldd[b1].Inflate(reader.hdr.GetNumberSamples(), TWK_LDD_LIST);
				blocks[0] = ldd[b1];
				blocks[0].Inflate(reader.hdr.GetNumberSamples(), unpack_list[method]);

				for(int i = 0; i < ldd[b1].n_rec; ++i){
					for(int j = i+1; j < ldd[b1].n_rec; ++j){
						if(blocks[0].blk->rcds[i].ac + blocks[0].blk->rcds[j].ac < 5){
							continue;
						}

						if(blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos > dist){
							//std::cerr << "out of range continue=" << (blocks[0].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) << std::endl;
							continue;
						}

						/*
						engine.PhasedVectorized(blocks[0], i, blocks[1], j);
						engine.PhasedVectorizedNoMissing(blocks[0], i, blocks[1], j);
						engine.PhasedVectorizedNoMissingNoTable(blocks[0], i, blocks[1], j);
						engine.UnphasedVectorized(blocks[0], i, blocks[1], j);
						engine.UnphasedVectorizedNoMissing(blocks[0], i, blocks[1], j);
						engine.Runlength(blocks[0], i, blocks[1], j);
						*/
						//engine.AllAlgorithms(blocks[0],i,blocks[0],j,&perfs[method]);
						(engine.*f[method])(blocks[0], i, blocks[0], j, &perfs[method]);
					}

				}
				n_total += ((blocks[0].n_rec * blocks[0].n_rec) - blocks[0].n_rec) / 2; // n choose 2

				for(int b2 = b1+1; b2 < reader.index.n; ++b2){
					//ldd[b2].Inflate(reader.hdr.GetNumberSamples(), TWK_LDD_LIST);
					blocks[1] = ldd[b2];
					blocks[1].Inflate(reader.hdr.GetNumberSamples(), unpack_list[method]);

					if(blocks[1].blk->rcds[0].pos - blocks[0].blk->rcds[blocks[0].n_rec-1].pos > dist){
						std::cerr << "never overlap=" << blocks[1].blk->rcds[0].pos - blocks[0].blk->rcds[blocks[0].n_rec-1].pos << std::endl;
						goto out_of_range;
					}

					for(int i = 0; i < blocks[0].n_rec; ++i){
						for(int j = 0; j < blocks[1].n_rec; ++j){
							if( blocks[0].blk->rcds[i].ac + blocks[1].blk->rcds[j].ac < 5 ){
								continue;
							}

							if(blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos > dist){
								//std::cerr << "out of range continue=" << (blocks[1].blk->rcds[j].pos - blocks[0].blk->rcds[i].pos) << std::endl;
								break;
							}
							/*
							engine.PhasedVectorized(blocks[0], i, blocks[1], j);
							engine.PhasedVectorizedNoMissing(blocks[0], i, blocks[1], j);
							engine.PhasedVectorizedNoMissingNoTable(blocks[0], i, blocks[1], j);
							engine.UnphasedVectorized(blocks[0], i, blocks[1], j);
							engine.UnphasedVectorizedNoMissing(blocks[0], i, blocks[1], j);
							engine.Runlength(blocks[0], i, blocks[1], j);
							*/
							(engine.*f[method])(blocks[0], i, blocks[1], j, &perfs[method]);
							//engine.AllAlgorithms(blocks[0],i,blocks[1],j,&perfs[method]);
						}
					}
					n_total += blocks[0].n_rec * blocks[1].n_rec;


					if(fl_lim++ == 10){
						std::cerr << "m=" << method << " " << utility::ToPrettyString(n_total) << " " << utility::ToPrettyString((uint64_t)((float)n_total/timer.Elapsed().count())) << "/s " << utility::ToPrettyString((uint64_t)((float)n_total*reader.hdr.GetNumberSamples()/timer.Elapsed().count())) << "/gt/s " << timer.ElapsedString() << " " << utility::ToPrettyString(engine.t_out) << std::endl;
						fl_lim = 0;
					}
				}
				out_of_range:
				continue;
			}

			std::cerr << "m=" << method << " " << utility::ToPrettyString(n_total) << " " << utility::ToPrettyString((uint64_t)((float)n_total/timer.Elapsed().count())) << "/s " << utility::ToPrettyString((uint64_t)((float)n_total*reader.hdr.GetNumberSamples()/timer.Elapsed().count())) << "/gt/s " << timer.ElapsedString() << " " << utility::ToPrettyString(engine.t_out) << std::endl;
		//}

		for(int i = 0; i < 10; ++i){
			for(int j = 0; j < perf_size; ++j){
				std::cout << i<<"\t"<<j<<"\t"<< (double)perfs[i].cycles[j]/perfs[i].freq[j] << "\t" << perfs[i].freq[j] << "\t" << perfs[i].cycles[j] << '\n';
			}
			delete[] perfs[i].cycles;
		}

		delete[] ldd;

		return true;
	}

public:
	twk_vimport_settings settings;

};

}



#endif /* IMPORTER_H_ */
