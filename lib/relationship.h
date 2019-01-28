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

void relationship_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Computes the relationship matrix for all samples.\n\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " relationship [options] -i <in.twk>\n\n"
	"Options:\n"
	"  -i FILE   input TWO file (required)\n"
	"  -I STRING interval string for target region\n\n";
}

int relationship(int argc, char** argv){
	if(argc < 3){
		relationship_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{"intervals",   required_argument, 0, 'I' },
		{"numeric",     optional_argument, 0, 'n' },
		{0,0,0,0}
	};

	tomahawk::twk_intervals ivals;
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
	std::cerr << tomahawk::utility::timestamp("LOG") << "Calling relationship..." << std::endl;

	tomahawk::twk_reader rdr;
	if(rdr.Open(input) == false){
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

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

	// Output matrix of size N^2.
	struct kinship_el {
		kinship_el() : n(0), c(0){}
		/*kinship_el& operator=(const kinship_el& other){
			n = other.n.load();
			c = other.c.load();
			return(*this);
		}*/

		//std::atomic<uint32_t> n, c;
		uint32_t n, c;
	};
	kinship_el** kin = new kinship_el*[rdr.hdr.GetNumberSamples()];
	for(int i = 0; i < rdr.hdr.GetNumberSamples(); ++i)
		kin[i] = new kinship_el[rdr.hdr.GetNumberSamples()];

	tomahawk::twk1_blk_iterator bit;
	bit.stream = rdr.stream;
	uint32_t n_variants = 0;

	// Todo: split matrix into N/t submatrices. One submatrix per worker thread.

	// If both are the same then add 1.0 -> 2
	// If both are different add 0.0 -> 0
	// If both have the same het layout at 0.5 -> 1
	// Divide end by 2
	tomahawk::Timer timer;
	for(int i = 0; i < n_blks; ++i){
		bit.stream->seekg(ivals.overlap_blocks[i]->foff);
		if(bit.NextBlock() == false){
			std::cerr << tomahawk::utility::timestamp("ERROR") << "Failed to load block " << i << "..." << std::endl;
			return false;
		}

		timer.Start();
		//n_variants += bit.blk.n;
		uint32_t score = 0;

		// Foreach record.
		for(int j = 0; j < bit.blk.n; ++j){
			if(ivals.itree[bit.blk.rcds[j].rid]->findOverlapping(bit.blk.rcds[j].pos, bit.blk.rcds[j].pos).size() == 0)
				continue;

			++n_variants;

			// Foreach RLE object against ever other (except self)
			// First loop over all, second loop over non-diagonal
			uint32_t sample_col = 0;
			for(int k = 0; k < bit.blk.rcds[j].gt->n; ++k){
				// Self
				//std::cerr << "outer=" << sample_col << "-" << sample_col + bit.blk.rcds[j].gt->GetLength(k) << std::endl;
				// Start a 1 because i==j then sample individual.
				for(int c = 0; c < bit.blk.rcds[j].gt->GetLength(k); ++c){
					for(int z = c + 1; z < bit.blk.rcds[j].gt->GetLength(k); ++z){
						// Add 2 for all these samples.
						kin[sample_col + c][sample_col + z].n += 2;
					}
				}

				uint32_t sample_row = sample_col + bit.blk.rcds[j].gt->GetLength(k);
				uint8_t refA = (bit.blk.rcds[j].gt->GetRefA(k) << 2) | bit.blk.rcds[j].gt->GetRefB(k);
				for(int l = k + 1; l < bit.blk.rcds[j].gt->n; ++l){
					uint8_t refB = (bit.blk.rcds[j].gt->GetRefA(l) << 2) | bit.blk.rcds[j].gt->GetRefB(l);
					// If both have different homozygous alleles then skip (add 0)
					if((refA == 0 && refB == 5) || (refA == 5 && refB == 0)){ // 0101 = 5
						// Add 0 for all these samples.
						sample_row += bit.blk.rcds[j].gt->GetLength(l);
						continue;
					}

					score = ((refA == 0 && refB == 0) || (refA == 5 && refB == 5) ? 2 : 1);

					// Sample col range (sample_col +c)
					for(int c = 0; c < bit.blk.rcds[j].gt->GetLength(k); ++c){
						// Sample row range (smaple_row + z)
						for(int z = 1; z < bit.blk.rcds[j].gt->GetLength(l); ++z){
							kin[sample_col + c][sample_row + z].n += score;
						}
					}
					sample_row += bit.blk.rcds[j].gt->GetLength(l);
				}
				assert(sample_row == rdr.hdr.GetNumberSamples());
				sample_col += bit.blk.rcds[j].gt->GetLength(k);
			}
			assert(sample_col == rdr.hdr.GetNumberSamples());
		}
		std::cerr << tomahawk::utility::timestamp("PROGRESS") << i << "/" << n_blks << "/" << rdr.index.n << ": " << bit.blk.n << " in " << timer.ElapsedString() << " (" << tomahawk::utility::ToPrettyString((uint64_t)((bit.blk.n*bit.blk.n*rdr.hdr.GetNumberSamples())/timer.Elapsed().count())) << " cmps/s)" << std::endl;
	}
	std::cerr << tomahawk::utility::timestamp("LOG") << "Number of individuals: " << rdr.hdr.GetNumberSamples() << " over " << n_variants << " sites..." << std::endl;

	// Add diagonal and add in lower triangular values.
	for(int i = 0; i < rdr.hdr.GetNumberSamples(); ++i) kin[i][i].n = 2*n_variants;
	for(int i = 0; i < rdr.hdr.GetNumberSamples(); ++i){
		for(int j = 1; j < rdr.hdr.GetNumberSamples(); ++j){
			kin[j][i] = kin[i][j];
		}
	}


	for(int i = 0; i < rdr.hdr.GetNumberSamples(); ++i){
		std::cout << (double)kin[i][0].n / n_variants / 2;
		for(int j = 1; j < rdr.hdr.GetNumberSamples(); ++j){
			std::cout << '\t' << (double)kin[i][j].n / n_variants / 2;
		}
		std::cout.put('\n');
	}
	std::cout.flush();

	for(int i = 0; i < rdr.hdr.GetNumberSamples(); ++i)
		delete[] kin[i];
	delete[] kin;

	return 0;
}
