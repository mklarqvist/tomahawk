#include "zstd.h"
#include "htslib/hts.h"
#include "tomahawk.h"

namespace tomahawk {

std::string LibrariesString(){
	return(std::string("Libraries: " + TOMAHAWK_PROGRAM_NAME + '-' + TOMAHAWK_LIB_VERSION
	              + "; ZSTD-" + std::string(ZSTD_versionString())
				  + "; htslib " + std::string(hts_version())));
}

void ProgramMessage(const bool separator){
	// General message for program version, git version, and linked library versions.
	std::cerr << "Program:   " << TOMAHAWK_PROGRAM_NAME << "-" << VERSION << " (Tools for computing, querying and storing LD data)" << std::endl;
	std::cerr << "Libraries: " << TOMAHAWK_PROGRAM_NAME << '-' << TOMAHAWK_LIB_VERSION
              << "; ZSTD-" << ZSTD_versionString()
			  << "; htslib " << std::string(hts_version()) << std::endl;
	std::cerr << "Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>" << std::endl;
	std::cerr << "Documentation: https://github.com/mklarqvist/tomahawk" << std::endl;
	std::cerr << "License: MIT" << std::endl;
	if(separator) std::cerr << "----------" << std::endl;
}

void ProgramHelp(void){
	std::cerr << "Usage: " << TOMAHAWK_PROGRAM_NAME << " [--version] [--help] <commands> <argument>" << std::endl;
	std::cerr << "Commands: aggregate, calc, scalc, concat, import, sort, view, haplotype, decay" << std::endl;
}

void ProgramHelpDetailed(void){
	ProgramHelp();
	std::cerr <<
    "\n"
	"aggregate    data rasterization framework for TWO files\n"
	"calc         calculate linkage disequilibrium\n"
	"scalc        calculate linkage disequilibrium for a single site\n"
	"concat       concatenate TWO files from the same set of samples\n"
	"import       import VCF/VCF.gz/BCF to TWK\n"
	"sort         sort TWO file\n"
    "view         TWO->LD/TWO view, TWO subset and filter\n"
	"haplotype    extract per-sample haplotype strings in FASTA/binary format\n"
	//"relationship compute marker-based pair-wise sample relationship matrices\n"
	"decay        compute LD-decay over distance\n"
	//"prune        perform graph-based LD-pruning of variant sites\n"
	//"stats        general stats for TWO files\n"
    << std::endl;
}

}
