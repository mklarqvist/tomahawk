#include "importer.h"

#include <thread>

#include "zstd.h"

#include "tomahawk.h"
#include "vcf_reader.h"
#include "buffer.h"
#include "genotype_encoder.h"
#include "core.h"
#include "index.h"
#include "zstd_codec.h"
#include "timer.h"
#include "writer.h"
#include "header_internal.h"

namespace tomahawk {

bool twk_variant_importer::Import(twk_vimport_settings& settings){
	this->settings = settings;
	return(this->Import());
}

bool twk_variant_importer::Import(void){
	// Start timer.
	Timer timer; timer.Start();

	if(settings.input != "-")
		std::cerr << utility::timestamp("LOG","READER") << "Opening " << settings.input << "..." << std::endl;

	// Retrieve a unique VcfReader.
	std::unique_ptr<VcfReader> vcf = tomahawk::VcfReader::FromFile(settings.input, std::thread::hardware_concurrency());
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
	if(settings.output.size() == 0 || (settings.output.size() == 1 && settings.output[0] == '-')){
		std::cerr << utility::timestamp("LOG","WRITER") << "Writing to stdout..." << std::endl;
		stream = &std::cout;
		stream_delete = false;
	}
	else {
		std::string base_path = tomahawk::twk_writer_t::GetBasePath(settings.output);
		std::string base_name = tomahawk::twk_writer_t::GetBaseName(settings.output);
		std::string extension = twk_writer_t::GetExtension(settings.output);
		if(extension.length() == 3){
			if(strncasecmp(&extension[0], "twk", 3) != 0){
				settings.output =  (base_path.size() ? base_path + "/" : "") + base_name + ".twk";
			}
		} else {
			 settings.output = (base_path.size() ? base_path + "/" : "") + base_name + ".twk";
		}

		std::cerr << utility::timestamp("LOG","WRITER") << "Opening " << settings.output << "..." << std::endl;
		stream = new std::ofstream;
		std::ofstream* outstream = reinterpret_cast<std::ofstream*>(stream);
		outstream->open(settings.output,std::ios::out | std::ios::binary);
		if(!outstream->good()){
			std::cerr << "failed to open" << std::endl;
			return false;
		}
	}

	// Append literal string.
	std::string import_string = "##tomahawk_importVersion=" + std::string(VERSION) + "\n";
	import_string += "##tomahawk_importCommand=" + tomahawk::LITERAL_COMMAND_LINE + "; Date=" + utility::datetime() + "\n";
	vcf->vcf_header_.literals_ += import_string;

	tomahawk::ZSTDCodec zcodec;
	stream->write(tomahawk::TOMAHAWK_MAGIC_HEADER.data(), tomahawk::TOMAHAWK_MAGIC_HEADER_LENGTH);
	tomahawk::twk_buffer_t buf(256000), obuf(256000);

	buf << vcf->vcf_header_;
	//std::cerr << "header buf size =" << buf.size() << std::endl;
	if(zcodec.Compress(buf, obuf, settings.c_level) == false){
		std::cerr << "failed to compress" << std::endl;
		return false;
	}
	//std::cerr << buf.size() << "->" << obuf.size() << " -> " << (float)buf.size()/obuf.size() << std::endl;

	stream->write(reinterpret_cast<const char*>(&buf.size()),sizeof(uint64_t));
	stream->write(reinterpret_cast<const char*>(&obuf.size()),sizeof(uint64_t));
	stream->write(obuf.data(),obuf.size());
	buf.reset();
	stream->flush();

	tomahawk::twk1_block_t block;
	uint32_t n_vnt_dropped = 0;
	uint32_t n_total_rec = 0;

	struct prev_helper {
		bool operator==(const bcf1_t* rec){
			if(rec->rid != rid) return false;
			if(rec->pos != pos) return false;
			if(dropped == true) return false;
			return true;
		}
		void operator=(const bcf1_t* rec){ dropped = false; rid = rec->rid; pos = rec->pos; }
		bool dropped = false;
		uint32_t rid, pos;
	};
	prev_helper prev_rec; prev_rec.rid = 0; prev_rec.pos = 0; prev_rec.dropped = false;

	memset(TWK_SITES_FILTERED, 0, sizeof(uint64_t)*8);
	uint64_t n_tot_vnts = 0;

	// While there are bcf1_t records available.
	while(vcf->next(BCF_UN_ALL)){
		++n_tot_vnts;
		tomahawk::twk1_t entry;
		entry.pos = vcf->bcf1_->pos;
		if(prev_rec == vcf->bcf1_){
			if(vcf->bcf1_->n_allele == 2){
				if(std::regex_match(vcf->bcf1_->d.allele[0],tomahawk::TWK_REGEX_CANONICAL_BASES) && std::regex_match(vcf->bcf1_->d.allele[1],tomahawk::TWK_REGEX_CANONICAL_BASES)){
					std::cerr << utility::timestamp("LOG") << "Duplicate site dropped: " << vcf->vcf_header_.GetContig(vcf->bcf1_->rid)->name << ":" << vcf->bcf1_->pos+1  << std::endl;
				}
			}
			prev_rec = vcf->bcf1_;
			prev_rec.dropped = true;
			++n_vnt_dropped;
			continue;
		}
		prev_rec = vcf->bcf1_;

		if(vcf->bcf1_->n_fmt){
			if(vcf->vcf_header_.GetFormat(vcf->bcf1_->d.fmt[0].id)->id == "GT"){
				if(vcf->bcf1_->d.fmt[0].n != 2){
					std::cerr << "do not support non-diploid samples" << std::endl;
					++TWK_SITES_FILTERED[3];
					++n_vnt_dropped;
					prev_rec.dropped = true;
					continue;
				}

				if(vcf->bcf1_->n_allele != 2){
					//std::cerr << "allele != 2" << std::endl;
					++TWK_SITES_FILTERED[6];
					++n_vnt_dropped;
					prev_rec.dropped = true;
					continue;
				}

				bool valid = true;
				for(int i = 0; i < 2; ++i){
					if(std::regex_match(vcf->bcf1_->d.allele[i],tomahawk::TWK_REGEX_CANONICAL_BASES) == false){
						//std::cerr << "failed match= " << vcf->bcf1_->d.allele[i] << std::endl;
						++TWK_SITES_FILTERED[7];
						++n_vnt_dropped;
						prev_rec.dropped = true;
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
				if(tomahawk::GenotypeEncoder::Encode(vcf->bcf1_, entry, settings) == false){
					//std::cerr << "invalid encoding" << std::endl;
					++n_vnt_dropped;
					prev_rec.dropped = true;
					continue;
				}
				entry.rid = vcf->bcf1_->rid;
				entry.calculateHardyWeinberg();

				if(entry.hwe < settings.hwe){
					//std::cerr << "dropped hwe=" << entry.hwe << "<" << settings.hwe << std::endl;
					++n_vnt_dropped;
					prev_rec.dropped = true;
					++TWK_SITES_FILTERED[8];
					continue;
				}

				if(block.n != 0){
					if(block.rid != vcf->bcf1_->rid){
						//std::cerr << "wrong rid@" << block.rid << "!=" << vcf->bcf1_->rid << std::endl;
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
						assert(zcodec.Compress(buf, obuf, settings.c_level));
						//std::cerr << buf.size() << "->" << obuf.size() << " = " << (float)buf.size()/obuf.size() << std::endl;
						// Write block: size-un, size, buffer
						tomahawk::twk_oblock_t oblock;
						oblock.Write(*stream, buf.size(), obuf.size(), obuf);

						ent.fend = stream->tellp();
						ent.b_unc = buf.size();
						ent.b_cmp = obuf.size();
						index += ent;
						//std::cerr << tomahawk::utility::timestamp("DEBUG") << outstream->tellp() << std::endl;
						buf.reset();
					}

					if(block.n == settings.block_size){
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
						ent.b_unc = buf.size();
						ent.b_cmp = obuf.size();
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
				prev_rec.dropped = true;
				++TWK_SITES_FILTERED[4];
				std::cerr << "no genotypes" << std::endl;
			}
		} else {
			++n_vnt_dropped;
			prev_rec.dropped = true;
			++TWK_SITES_FILTERED[5];
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
		assert(zcodec.Compress(buf, obuf, settings.c_level));
		//std::cerr << buf.size() << "->" << obuf.size() << " = " << (float)buf.size()/obuf.size() << std::endl;
		// Write block: size-un, size, buffer
		tomahawk::twk_oblock_t oblock;
		oblock.Write(*stream, buf.size(), obuf.size(), obuf);

		ent.fend = stream->tellp();
		ent.b_unc = buf.size();
		ent.b_cmp = obuf.size();
		index += ent;

		//std::cerr << tomahawk::utility::timestamp("DEBUG") << "Final=" << outstream->tellp() << std::endl;
		buf.reset();
	}

	buf << index;
	//std::cerr << "index buf size =" << buf.size() << std::endl;
	if(zcodec.Compress(buf, obuf, settings.c_level) == false){
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
	std::cerr << utility::timestamp("LOG") << "Filtered out " << utility::ToPrettyString(n_tot_vnts - index.GetTotalVariants()) << " sites (" << (float)(n_tot_vnts - index.GetTotalVariants())/n_tot_vnts*100 << "%):" << std::endl;
	for(int i = 0; i < 9; ++i){
		std::cerr << utility::timestamp("LOG") << "   " << TWK_SITES_FILTERED_NAMES[i] << ": " << utility::ToPrettyString(TWK_SITES_FILTERED[i]) << " (" << (float)TWK_SITES_FILTERED[i]/n_tot_vnts*100 << "%)" << std::endl;
	}

	if(stream_delete) delete stream;
	return(true);
}
}
