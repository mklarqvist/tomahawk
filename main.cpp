#include <fstream>

#include "tomahawk.h"
std::string tomahawk::LITERAL_COMMAND_LINE;
std::string tomahawk::INTERPRETED_COMMAND;

#include "importer.h"
#include "ld.h"
#include "two_reader.h"

int main(void){
	if(0){
		tomahawk::twk_vimport_settings settings;
		tomahawk::twk_variant_importer importer;
		settings.input  = "/home/mk21/Downloads/1kgp3_chr20.bcf";
		settings.output = "/home/mk21/Downloads/test.twk";
		if(importer.Import(settings) == false){
			std::cerr << "failed import" << std::endl;
			return 1;
		}
	}

	if(1){
		tomahawk::twk_ld_settings settings;
		tomahawk::twk_ld ld;
		settings.in  = "/home/mk21/Downloads/test.twk";
		settings.out = "/home/mk21/Downloads/debug.two";
		settings.n_chunks = 1; // has to be non-negative, non-zero
		settings.c_chunk  = 0; // has to be in range [0, n_chunks)
		settings.minR2    = 0.1;
		settings.minP     = 1e-3;
		settings.window   = true;
		settings.l_window = 500000;

		if(ld.Compute(settings) == false){
			std::cerr << "failed compute" << std::endl;
			return false;
		}
	}

	tomahawk::two_reader oreader;
	if(oreader.Open("/home/mk21/Downloads/debug.two") == false){
		std::cerr << "failed to open" << std::endl;
		return false;
	}

	tomahawk::twk1_two_iterator it;
	it.stream = oreader.stream;
	uint32_t i = 0;
	tomahawk::twk1_two_block_t blk2;
	std::cerr << "resizing to=" << 1000000000/sizeof(tomahawk::twk1_two_t) << " entries..." << std::endl;
	blk2.resize(1000000000/sizeof(tomahawk::twk1_two_t));

	uint32_t tot = 0;
	while(it.NextBlock()){
		//std::cerr << i++ << " -> " << it.GetBlock().n << std::endl;
		tot += it.GetBlock().n;
		for(int j = 0; j < it.blk.n; ++j){
			//it.blk[j].Print(std::cerr);
			if(blk2.n == blk2.m){
				//std::cerr << "resseting" << std::endl;
				blk2.Sort();
				//for(int k = 0; k < blk2.n; ++k){
				//	blk2.rcds[k].Print(std::cerr);
				//}
				blk2.reset();
			}
			blk2 += it.blk[j];
			//it.blk[j].Print(std::cout);
		}
		//std::cerr << blk2.n << "/" << blk2.m << std::endl;
	}
	if(blk2.n){
		//std::cerr << "resseting" << std::endl;
		blk2.Sort();
		//for(int k = 0; k < blk2.n; ++k){
		//	blk2.rcds[k].Print(std::cerr);
		//}
		blk2.reset();
	}
	std::cerr << "total=" << tot << "/" << oreader.index.GetTotalVariants() << std::endl;

	return 0;
}
