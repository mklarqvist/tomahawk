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
		settings.input  = "/media/mdrk/NVMe/1kgp3/bcf/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf";
		settings.output = "/media/mdrk/NVMe/1kgp3/test.twk";
		if(importer.Import(settings) == false){
			std::cerr << "failed import" << std::endl;
			return 1;
		}

		tomahawk::twk_ld ld;
		ld.Compute();
	}

	tomahawk::two_reader oreader;
	if(oreader.Open("/media/mdrk/NVMe/1kgp3/debug.two") == false){
		std::cerr << "failed to open" << std::endl;
		return false;
	}
	tomahawk::twk1_two_iterator it;
	it.stream = oreader.stream;
	uint32_t i = 0;
	tomahawk::twk1_two_block_t blk2;
	blk2.resize(1000000);

	uint32_t tot = 0;
	while(it.NextBlock()){
		//std::cerr << i++ << " -> " << it.GetBlock().n << std::endl;
		tot += it.GetBlock().n;
		for(int j = 0; j < it.blk.n; ++j){
			//it.blk[j].Print(std::cerr);
			if(blk2.n == blk2.m){
				std::cerr << "resseting" << std::endl;
				blk2.Sort();
				//for(int k = 0; k < blk2.n; ++k){
				//	blk2.rcds[k].Print(std::cerr);
				//}
				blk2.reset();
			}
			blk2 += it.blk[j];
		}
		//std::cerr << blk2.n << "/" << blk2.m << std::endl;
	}
	std::cerr << "total=" << tot << std::endl;

	return 0;
}
