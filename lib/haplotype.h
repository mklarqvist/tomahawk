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
#include "twk_reader.h"

void haplotype_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Retrieve haplotype strings from individuals in a given range.\n"
	"        This is equivalent to slicing out a section of the transpose of the\n"
	"        genotype-sample matrix. Outputs a FASTA file with 2N haplotype strings.\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " haplotype [options] -i <in.twk> -I <interval>\n\n"
	"Options:\n"
	"  -i FILE   input TWO file (required)\n"
	"  -I STRING interval string for target region\n\n";
}

int haplotype(int argc, char** argv){
	if(argc < 3){
		haplotype_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{"intervals",   required_argument, 0, 'I' },
		{"numeric",   optional_argument, 0, 'n' },
		{0,0,0,0}
	};

	tomahawk::twk_intervals_twk ivals;
	std::string input;
	std::vector<std::string> intervals;
	bool output_numeric_encoding = false;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:I:n?", long_options, &long_index)) != -1){
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
			input = std::string(optarg);
			break;
		case 'I':
			intervals.push_back(std::string(optarg));
			break;
		case 'n':
			output_numeric_encoding = true;
			break;
		}
	}

	if(input.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(intervals.size() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No interval(s) provided..." << std::endl;
		return(1);
	}

	// Print messages
	tomahawk::ProgramMessage();
	std::cerr << tomahawk::utility::timestamp("LOG") << "Calling haplotype..." << std::endl;

	tomahawk::twk_reader rdr;
	if(rdr.Open(input)){
		ivals.ivecs.resize(rdr.hdr.GetNumberContigs());
		if(ivals.ParseIntervalStrings(intervals, rdr.hdr) == false){
			std::cerr << "failed to parse intervals" << std::endl;
			return 1;
		}

		if(ivals.Build(rdr.hdr.GetNumberContigs(), rdr.index) == false){
			std::cerr << "failed to build" << std::endl;
			return 1;
		}

		// Number of blocks available to iterate over.
		const uint32_t n_blks = ivals.overlap_blocks.size();
		if(n_blks == 0){
			std::cerr << "no data available in that range" << std::endl;
			return 1;
		}

		// Output matrix.
		std::vector< std::vector<uint8_t> > haps(2*rdr.hdr.GetNumberSamples(), std::vector<uint8_t>());

		tomahawk::twk1_blk_iterator bit;
		bit.stream = rdr.stream;
		uint32_t n_variants = 0;
		char fasta_lookup[3];
		fasta_lookup[2] = 'N';
		char numeric_lookup[3];
		numeric_lookup[0] = '0', numeric_lookup[1] = '1', numeric_lookup[2] = '2';
		char* lookup = fasta_lookup;
		if(output_numeric_encoding){
			std::cerr << "using numeric lookup" << std::endl;
			lookup = numeric_lookup;
		}

		for(int i = 0; i < n_blks; ++i){
			bit.stream->seekg(ivals.overlap_blocks[i]->foff);
			if(bit.NextBlock() == false){
				std::cerr << tomahawk::utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
				return false;
			}

			std::cerr << i << "/" << n_blks << "/" << rdr.index.n << ": " << bit.blk.n << std::endl;
			n_variants += bit.blk.n;
			// Foreach record.
			for(int j = 0; j < bit.blk.n; ++j){
				//std::cerr << bit.blk.rcds[j].pos << "\t" << bit.blk.rcds[j].gt->GetLength(0);
				// Foreach RLE object.
				uint32_t hap_offset = 0;
				fasta_lookup[0] = bit.blk.rcds[j].GetAlleleA();
				fasta_lookup[1] = bit.blk.rcds[j].GetAlleleB();
				for(int k = 0; k < bit.blk.rcds[j].gt->n; ++k){
					for(int p = 0; p < bit.blk.rcds[j].gt->GetLength(k); ++p, hap_offset += 2){
						//std::cerr << "," << bit.blk.rcds[j].gt->GetLength(k);
						haps[hap_offset].push_back(lookup[bit.blk.rcds[j].gt->GetRefA(p)]);
						haps[hap_offset+1].push_back(lookup[bit.blk.rcds[j].gt->GetRefB(p)]);
					}
				}
				//std::cerr << std::endl;
			}
		}
		std::cerr << haps.size() << "," << haps[0].size() << std::endl;
		for(int p = 0; p < haps.size(); ++p){
			std::cout << ">" << rdr.hdr.samples_[p/2] << "_" << (p%2) << "\n";
			for(int i = 0; i < haps[p].size(); ++i){
				std::cout << haps[p][i];
			}
			std::cout.put('\n');
		}
		std::cout.flush();

		std::cerr << "done. variants=" << n_variants << std::endl;
		return 0;
	}
	else {
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

	//tomahawk::two_sorter sorter;
	//if(sorter.Sort(settings) == false) return 1;
	//return 0;
}
