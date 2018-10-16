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
		// Retrieve a unique VcfReader.
		std::string filename = "/home/mk21/Downloads/1kgp3_chr20.bcf";
		//filename = "/home/mk21/Downloads/fish_callset_filtered.bcf";
		filename = "/home/mk21/Downloads/randomized.100k-50k.bcf";
		std::unique_ptr<tomahawk::io::VcfReader> vcf = tomahawk::io::VcfReader::FromFile(filename, 8);
		if(vcf == nullptr){
			std::cerr << "failed to get vcfreader" << std::endl;
			return false;
		}
		if(vcf->vcf_header_.GetFormat("GT") == nullptr){
			std::cerr << "Genotype data not set in this file" << std::endl;
			return false;
		}
		// Todo: add header literal tracing input parameters

		tomahawk::Index index(vcf->vcf_header_.GetNumberContigs());

		std::ostream* stream = new std::ofstream;
		std::ofstream* outstream = reinterpret_cast<std::ofstream*>(stream);
		outstream->open("/home/mk21/Downloads/debug.twk",std::ios::out | std::ios::binary);
		if(!outstream->good()){
			std::cerr << "failed to open" << std::endl;
			return false;
		}

		tomahawk::ZSTDCodec zcodec;
		outstream->write(tomahawk::TOMAHAWK_MAGIC_HEADER.data(), tomahawk::TOMAHAWK_MAGIC_HEADER_LENGTH);
		tomahawk::twk_buffer_t buf(256000), obuf(256000);
		buf << vcf->vcf_header_;
		std::cerr << "header buf size =" << buf.size() << std::endl;
		assert(zcodec.Compress(buf, obuf, 10));
		std::cerr << buf.size() << "->" << obuf.size() << " -> " << (float)buf.size()/obuf.size() << std::endl;

		outstream->write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		outstream->write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		outstream->write(obuf.data(),obuf.size());
		buf.reset();

		// start
		outstream->flush();
		const uint64_t start_pos = outstream->tellp();

		tomahawk::twk1_block_t block;
		uint32_t n_vnt_dropped = 0;


		while(vcf->next(BCF_UN_ALL)){
			tomahawk::twk1_t entry;
			entry.pos = vcf->bcf1_->pos;

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
							std::cerr << tomahawk::utility::timestamp("DEBUG") << outstream->tellp() << std::endl;
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
							std::cerr << buf.size() << "->" << obuf.size() << " = " << (float)buf.size()/obuf.size() << std::endl;
							// Write block: size-un, size, buffer
							tomahawk::twk_oblock_t oblock;
							oblock.Write(*stream, buf.size(), obuf.size(), obuf);

							ent.fend = stream->tellp();
							index += ent;

							std::cerr << tomahawk::utility::timestamp("DEBUG") << outstream->tellp() << std::endl;
							buf.reset();
						}
						block += entry;
					} else {
						block += entry;
						block.rid = vcf->bcf1_->rid;
						block.minpos = vcf->bcf1_->pos;
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
			std::cerr << buf.size() << "->" << obuf.size() << " = " << (float)buf.size()/obuf.size() << std::endl;
			// Write block: size-un, size, buffer
			tomahawk::twk_oblock_t oblock;
			oblock.Write(*stream, buf.size(), obuf.size(), obuf);

			ent.fend = stream->tellp();
			index += ent;

			std::cerr << tomahawk::utility::timestamp("DEBUG") << "Final=" << outstream->tellp() << std::endl;
			buf.reset();
		}

		buf << index;
		std::cerr << "index buf size =" << buf.size() << std::endl;
		assert(zcodec.Compress(buf, obuf, 10));
		std::cerr << buf.size() << "->" << obuf.size() << " -> " << (float)buf.size()/obuf.size() << std::endl;
		const uint64_t offset_start_index = stream->tellp();
		uint8_t marker = 0;
		outstream->write(reinterpret_cast<const char*>(&marker),sizeof(uint8_t));
		outstream->write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
		outstream->write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
		outstream->write(obuf.data(),obuf.size());
		outstream->write(reinterpret_cast<const char*>(&offset_start_index),sizeof(uint64_t));
		outstream->write(tomahawk::TOMAHAWK_FILE_EOF.data(), tomahawk::TOMAHAWK_FILE_EOF_LENGTH);
		outstream->flush();

		delete stream;

		//*//////////////// Repoen and compute
		twk_reader reader;
		if(reader.Open("/home/mk21/Downloads/debug.twk") == false){
			std::cerr << "failed" << std::endl;
			return false;
		}

		twk1_blk_iterator bit;
		bit.stream = reader.rstream;
		LDEngine engine; engine.SetSamples(vcf->vcf_header_.GetNumberSamples());

		std::cerr << "index entries=" << index.n << std::endl;

		Timer timer; timer.Start();
		twk1_ldd_blk* ldd = new twk1_ldd_blk[index.n];

		uint32_t i = 0;
		while(bit.NextBlockRaw()){
			std::cerr << "next-block=" << bit.stream->tellg() << "/" << offset_start_index << " n=" << bit.blk.n << std::endl;
			ldd[i++].SetOwn(bit, vcf->vcf_header_.GetNumberSamples());
		}
		std::cerr << "Done = " << timer.ElapsedString() << std::endl;
		timer.Start();

		uint64_t n_total = 0; uint32_t fl_lim = 0;
		for(int b1 = 0; b1 < index.n; ++b1){
			ldd[b1].Inflate(vcf->vcf_header_.GetNumberSamples());
			for(int i = 0; i < ldd[b1].n_rec; ++i){
				for(int j = i+1; j < ldd[b1].n_rec; ++j){
					if( ldd[b1].blk->rcds[i].ac + ldd[b1].blk->rcds[j].ac < 5){
						continue;
					}
					/*
					if( std::max(ldd[b1].blk->rcds[i].ac, ldd[b1].blk->rcds[j].ac) > 1000 ){
						engine.PhasedVectorizedNoMissingNoTable(ldd[b1].vec[i], ldd[b1].vec[j]);
					}
					else
					*/
					engine.CalculateLDPhasedSimple(ldd[b1].list[i], ldd[b1].list[j]);

				}

			}
			n_total += ((ldd[b1].n_rec * ldd[b1].n_rec) - ldd[b1].n_rec) / 2; // n choose 2

			for(int b2 = b1+1; b2 < index.n; ++b2){
				ldd[b2].Inflate(vcf->vcf_header_.GetNumberSamples());
				for(int i = 0; i < ldd[b1].n_rec; ++i){
					for(int j = 0; j < ldd[b2].n_rec; ++j){
						if( ldd[b1].blk->rcds[i].ac + ldd[b1].blk->rcds[j].ac < 5){
							continue;
						}
						/*
						if( std::max(ldd[b1].blk->rcds[i].ac, ldd[b1].blk->rcds[j].ac) > 1000 ){
							engine.PhasedVectorizedNoMissingNoTable(ldd[b1].vec[i], ldd[b1].vec[j]);
						}
						else
						*/
						engine.CalculateLDPhasedSimple(ldd[b1].list[i], ldd[b1].list[j]);
					}
				}
				n_total += ldd[b1].n_rec * ldd[b2].n_rec;
				if(fl_lim++ == 10){
					std::cerr << "DoneLD = " << utility::ToPrettyString(n_total) << " " << utility::ToPrettyString((uint64_t)((float)n_total/timer.Elapsed().count())) << "/s " << timer.ElapsedString() << std::endl;
					fl_lim = 0;
				}
				ldd[b2].Clear();
			}
			ldd[b1].Clear();
		}

		std::cerr << "DoneLD = " << utility::ToPrettyString(n_total) << " " << utility::ToPrettyString((uint64_t)((float)n_total/timer.Elapsed().count())) << "/s " << timer.ElapsedString() << std::endl;


		//std::cerr << "done=" << bit.blk.n*bit2.blk.n << std::endl;

		delete[] ldd;

		return true;
	}

public:
	twk_vimport_settings settings;

};

}



#endif /* IMPORTER_H_ */
