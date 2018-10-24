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
		//if(importer.Import(settings) == false){
		//	std::cerr << "failed import" << std::endl;
		//	return 1;
		//}

		tomahawk::twk_ld ld;
		ld.Compute();
	}

	tomahawk::two_reader oreader;
	if(oreader.Open("/home/mk21/Downloads/debug.two") == false){
		std::cerr << "failed to open" << std::endl;
		return false;
	}
	tomahawk::twk1_two_iterator it;
	it.stream = oreader.stream;
	uint32_t i = 0;
	while(it.NextBlock()){
		std::cerr << i++ << " -> " << it.GetBlock().n << std::endl;
		for(int j = 0; j < it.blk.n; ++j){
			it.blk[j].Print(std::cerr);
		}
	}

	return 0;
}
