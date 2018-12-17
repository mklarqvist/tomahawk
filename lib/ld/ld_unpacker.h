#ifndef LIB_LD_LD_UNPACKER_H_
#define LIB_LD_LD_UNPACKER_H_

#include "core.h"
#include "twk_reader.h"

#include <thread>

namespace tomahawk {

/**<
 * Supportive structure that unpacks RLE-encoded genotypes into the appropriate
 * data structure for use in calculating LD. Most commonly, RLE vectors are
 * unpacked into either bitvectors or bitmaps.
 */
struct twk_ld_unpacker {
	/**<
	 * Start unpacking by spawning a std::thread performing the computation
	 * over the ranges provided. Spawned threads has to be joined _outside_
	 * this function. No checks are made.
	 * @param diag     Is the data diagonal?
	 * @param settings Reference instance of the user-settings object.
	 * @return         Returns a pointer to the spawned std::thread or nullptr if failed.
	 */
	std::thread* Start(const bool diag, const int32_t load_type, const std::string& in){
		load = load_type;
		this->diag = diag;

		bit.stream = new std::ifstream(in, std::ios::binary | std::ios::in);
		if(bit.stream->good() == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to open \"" << in << "\"!" << std::endl;
			return nullptr;
		}

		if(diag) thread = new std::thread(&twk_ld_unpacker::UnpackDiagonal, this);
		else thread = new std::thread(&twk_ld_unpacker::UnpackSquare, this);
		return(thread);
	}

	/**<
	 * Unpacks data from a diagonal block.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool UnpackDiagonal(){
		bit.stream->seekg(rdr->index.ent[fL].foff);
		if(bit.stream->good() == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to seek to index offset " << fL << " -> " << rdr->index.ent[fL].foff << "!" << std::endl;
			return false;
		}

		//std::cerr << "in unpack diag (" << fL << "-" << tL << ")" << std::endl;
		//std::cerr << "left-shift=" << lshift << std::endl;
		for(int i = lshift; i < lshift + (tL-fL); ++i){
			//std::cerr << "at1=" << i << "/" << lshift + (tL-fL) << std::endl;
			if(bit.NextBlock() == false){
				std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "!" << std::endl;
				return false;
			}

			ldd2[i] = std::move(bit.blk);
			ldd[i].SetOwn(ldd2[i], rdr->hdr.GetNumberSamples());
			ldd[i].Inflate(rdr->hdr.GetNumberSamples(),load, resize);
		}

		std::ifstream* s = reinterpret_cast<std::ifstream*>(bit.stream);
		s->close();
		bit.stream = nullptr;

		return true;
	}

	/**<
	 * Unpack data from square blocks (blocks never overlap).
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool UnpackSquare(){
		bit.stream->seekg(rdr->index.ent[fL].foff);
		if(bit.stream->good() == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to seek to index offset " << fL << " -> " << rdr->index.ent[fL].foff << "!" << std::endl;
			return false;
		}

		//std::cerr << "in unpack square (" << fL << "-" << tL << ")(" << fR << "-" << tR << ")" << std::endl;
		//std::cerr << "left-shift=" << lshift << std::endl;
		for(int i = lshift; i < lshift + (tL-fL); ++i){
			//std::cerr << "at1=" << i << "/" << lshift + (tL-fL) << std::endl;
			if(bit.NextBlock() == false){
				std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "!" << std::endl;
				return false;
			}

			ldd2[i] = std::move(bit.blk);
			ldd[i].SetOwn(ldd2[i], rdr->hdr.GetNumberSamples());
			ldd[i].Inflate(rdr->hdr.GetNumberSamples(),load, resize);
		}

		bit.stream->seekg(rdr->index.ent[fR].foff); // seek absolute offset
		if(bit.stream->good() == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to seek to index offset " << fR << " -> " << rdr->index.ent[fR].foff << "!" << std::endl;
			return false;
		}

		// compute with local offset
		//std::cerr << "left offset=" << loff << " and right=" << roff << std::endl;
		//std::cerr << "offset=" << fR - loff << " and steps=" << loff + (tR - fR) << std::endl;
		for(int i = loff + roff; i < loff + roff + (tR - fR); ++i){
			//std::cerr << "at2=" << i << "/" << loff + roff + (tR - fR) << " with range=" << (tR - fR) << std::endl;
			if(bit.NextBlock() == false){
				std::cerr << utility::timestamp("ERROR") << "Failed to load block " << i << "!" << std::endl;
				return false;
			}

			ldd2[i] = std::move(bit.blk);
			ldd[i].SetOwn(ldd2[i], rdr->hdr.GetNumberSamples());
			ldd[i].Inflate(rdr->hdr.GetNumberSamples(),load, resize);
		}

		std::ifstream* s = reinterpret_cast<std::ifstream*>(bit.stream);
		s->close();
		bit.stream = nullptr;

		return true;
	}

public:
	bool resize, diag;
	uint8_t load;
	uint32_t fL, tL, fR, tR, loff, roff, lshift;
	twk_reader* rdr;
	std::thread* thread;
	twk1_ldd_blk* ldd;
	twk1_block_t* ldd2;
	twk1_blk_iterator bit;
};

}



#endif /* LIB_LD_LD_UNPACKER_H_ */
