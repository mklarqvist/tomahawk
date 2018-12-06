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
#include "two_reader.h"

struct sstats {
	typedef void (sstats::*addfunc)(const tomahawk::twk1_two_t*);
	typedef double (sstats::*redfunc)(const uint32_t) const;

	sstats() : n(0), total(0), total_squared(0), min(0), max(0){}

	void AddR2(const tomahawk::twk1_two_t* rec){ Add(rec->R2); }
	void AddR(const tomahawk::twk1_two_t* rec){ Add(rec->R); }
	void AddD(const tomahawk::twk1_two_t* rec){ Add(rec->D); }
	void AddDprime(const tomahawk::twk1_two_t* rec){ Add(rec->Dprime); }
	void AddP(const tomahawk::twk1_two_t* rec){ Add(rec->P); }
	void AddHets(const tomahawk::twk1_two_t* rec){
		Add((rec->cnt[1] + rec->cnt[2]) / (rec->cnt[0] + rec->cnt[1] + rec->cnt[2] + rec->cnt[3]));
	}

	void AddAlts(const tomahawk::twk1_two_t* rec){
		Add((rec->cnt[3]) / (rec->cnt[0] + rec->cnt[1] + rec->cnt[2] + rec->cnt[3]));
	}

	double GetMean(const uint32_t min = 0) const {
		if(n < min) return(0);
		return(total / n);
	}

	double GetCount(const uint32_t min = 0) const {
		if(n < min) return(0);
		return(n);
	}

	template <class T> void Add(const T value, const double weight = 1){
		this->total         += value;
		this->total_squared += value*value;
		this->n             += weight;
		if(value < this->min) this->min = value;
		if(value > this->max) this->max = value;
	}

	double GetStandardDeviation(const uint32_t cutoff = 2) const{
		if(this->n < cutoff) return(0);
		return(sqrt(this->total_squared/this->n - (this->total / this->n)*(this->total / this->n)));
	}

	// Accessor functions
	inline double GetTotal(const uint32_t cutoff = 0) const{ return(this->total); }
	inline double GetTotalSquared(const uint32_t cutoff = 0) const{ return(this->total_squared); }
	inline double GetMin(const uint32_t cutoff = 0) const{ return(this->min); }
	inline double GetMax(const uint32_t cutoff = 0) const{ return(this->max); }

	uint64_t n;
	double total, total_squared;
	double min, max;
};

void aggregate_usage(void){
	tomahawk::ProgramMessage();
	std::cerr <<
	"About:  Aggregate TWO data into a rasterized matrix of size [x,y] for\n"
    "        plotting. Data has to be sorted and pre-sliced to the correct\n"
	"        interval-of-interest prior to running.\n"
	"Usage:  " << tomahawk::TOMAHAWK_PROGRAM_NAME << " aggregate [options] <in.two>\n\n"

	"Options:\n"
	"  -i   FILE   input TWO file (required)\n"
	"  -x,y INT    number of X/Y-axis bins (default: 1000)\n"
	"  -f   STRING aggregation function: can be one of (r2,r,d,dprime,dp,p,hets,alts,het,alt)(required)\n"
	"  -r   STRING reduction function: can be one of (mean,count,n,min,max,sd)(required)\n"
	"  -I   STRING filter interval <contig>:pos-pos (TWK/TWO) or linked interval <contig>:pos-pos,<contig>:pos-pos\n"
	"  -c   INT    min cut-off value used in reduction function: value < c will be set to 0 (default: 5)\n\n";
}

int aggregate(int argc, char** argv){
	if(argc < 3){
		aggregate_usage();
		return(0);
	}

	static struct option long_options[] = {
		{"input",       required_argument, 0, 'i' },
		{"interval",    optional_argument, 0, 'I' },
		{"xbins",    optional_argument, 0, 'x' },
		{"ynins",    optional_argument, 0, 'y' },
		{"min-cutoff",    optional_argument, 0, 'c' },
		{"aggregate-function",    required_argument, 0, 'f' },
		{"reduce-function",    required_argument, 0, 'r' },

		{0,0,0,0}
	};

	tomahawk::twk_two_settings settings;
	int32_t x_bins = 1000, y_bins = 1000;
	std::string aggregate_func_name, reduce_func_name;
	int32_t min_cutoff = 5;

	int c = 0;
	int long_index = 0;
	int hits = 0;
	while ((c = getopt_long(argc, argv, "i:I:x:y:f:r:c:", long_options, &long_index)) != -1){
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
			settings.in = std::string(optarg);
			break;
		case 'I': settings.ivals.push_back(std::string(optarg)); break;
		case 'x': x_bins = std::atoi(optarg); break;
		case 'y': y_bins = std::atoi(optarg); break;
		case 'f': aggregate_func_name = std::string(optarg); break;
		case 'r': reduce_func_name = std::string(optarg); break;
		case 'c': min_cutoff = std::atoi(optarg); break;
		}
	}

	sstats::addfunc f = &sstats::AddR2;
	sstats::redfunc r = &sstats::GetMean;

	if(aggregate_func_name.size() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No aggregation function (-f) provided..." << std::endl;
		return(1);
	}

	if(reduce_func_name.size() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No reduce function (-r) provided..." << std::endl;
		return(1);
	}

	if(min_cutoff < 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Cannot have a min-cutoff (-c) < 0..." << std::endl;
		return(1);
	}

	std::transform(aggregate_func_name.begin(), aggregate_func_name.end(), aggregate_func_name.begin(), ::tolower);
	if(aggregate_func_name == "r2"){ f = &sstats::AddR2; }
	else if(aggregate_func_name == "r"){ f = &sstats::AddR; }
	else if(aggregate_func_name == "d"){ f = &sstats::AddD; }
	else if(aggregate_func_name == "dprime"){ f = &sstats::AddDprime; }
	else if(aggregate_func_name == "dp"){ f = &sstats::AddDprime; }
	else if(aggregate_func_name == "p"){ f = &sstats::AddP; }
	else if(aggregate_func_name == "hets"){ f = &sstats::AddHets; }
	else if(aggregate_func_name == "alts"){ f = &sstats::AddAlts; }
	else if(aggregate_func_name == "het"){ f = &sstats::AddHets; }
	else if(aggregate_func_name == "alt"){ f = &sstats::AddAlts; }
	else {
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Unknown aggregation function \"" << aggregate_func_name << "\"..." << std::endl;
		return(1);
	}

	std::transform(reduce_func_name.begin(), reduce_func_name.end(), reduce_func_name.begin(), ::tolower);
	if(reduce_func_name == "mean"){ r = &sstats::GetMean; }
		else if(reduce_func_name == "max"){ r = &sstats::GetMax; }
		else if(reduce_func_name == "min"){ r = &sstats::GetMin; }
		else if(reduce_func_name == "count"){ r = &sstats::GetTotal; }
		else if(reduce_func_name == "n"){ r = &sstats::GetTotal; }
		else if(reduce_func_name == "sd"){ r = &sstats::GetStandardDeviation; }
		else {
			std::cerr << tomahawk::utility::timestamp("ERROR") << "Unknown reduce function \"" << reduce_func_name << "\"..." << std::endl;
			return(1);
		}

	if(x_bins < 5){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Number of x-bins cannot be < 5" << std::endl;
		return(1);
	}

	if(y_bins < 5){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "Number of y-bins cannot be < 5" << std::endl;
		return(1);
	}

	if(settings.in.length() == 0){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	// New instance of reader.
	tomahawk::two_reader oreader;

	// Open file handle.
	if(oreader.Open(settings.in) == false){
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

	// Build intervals data structures if any are available.
	if(settings.intervals.Build(settings.ivals,
	                            oreader.hdr.GetNumberContigs(),
	                            oreader.index,
	                            oreader.hdr) == false)
	{
		return 1;
	}

	// Construct filters.
	settings.filter.Build();

	// Algorithmic overview.
	// Step 1: Find maximum and minimum X and Y values.
	// Step 2: Partition into (maxY-minY)/#bins and (maxX-minX)/#bins buckets.
	//      a: Range can be either dynamic min and max given the data, or;
	//      b: Interval (from,to)-tuple.
	// Step 3: Iterate over data and update summary statistics in buckets.
	// Step 4: Output data.
	if(oreader.index.state != TWK_IDX_SORTED){
		std::cerr << tomahawk::utility::timestamp("ERROR") << "The input file has to be sorted..." << std::endl;
		return 1;
	}

	// Step 1: Calculate the landscape ranges (X and Y dimensions).
	// Approach 1: Without dropping regions with no data.
	uint64_t range = 0;
	for(int i = 0; i < oreader.index.m_ent; ++i){
		std::cerr << oreader.index.ent_meta[i].rid << ":" << oreader.index.ent_meta[i].minpos << "-" << oreader.index.ent_meta[i].maxpos << " and n=" << oreader.index.ent_meta[i].n << "," << oreader.index.ent_meta[i].nn << std::endl;
		if(oreader.index.ent_meta[i].n)
			range += (oreader.index.ent_meta[i].maxpos - oreader.index.ent_meta[i].minpos) + 1;
	}
	std::cerr << "range=" << range << std::endl;

	// Approach 2: Dropping regions with no data.
	struct offset_tuple {
		offset_tuple() : range(0), min(0), max(0){}
		uint64_t range;
		uint32_t min, max;
	};
	std::vector<offset_tuple> rid_offsets(oreader.index.m_ent);
	for(int i = 0; i < oreader.index.m_ent; ++i){
		std::cerr << oreader.index.ent_meta[i].rid << ":" << oreader.index.ent_meta[i].minpos << "-" << oreader.index.ent_meta[i].maxpos << " and n=" << oreader.index.ent_meta[i].n << "," << oreader.index.ent_meta[i].nn << std::endl;
		if(oreader.index.ent_meta[i].n){
			// Cumulative offset for the current rid equals the previous rid
			if(oreader.index.ent_meta[i].rid != 0){
				rid_offsets[oreader.index.ent_meta[i].rid].range = rid_offsets[oreader.index.ent_meta[i].rid - 1].range + ((oreader.index.ent_meta[i].maxpos - oreader.index.ent_meta[i].minpos) + 1);
				range += (oreader.index.ent_meta[i].maxpos - oreader.index.ent_meta[i].minpos) + 1;
			}

		} else {
			if(oreader.index.ent_meta[i].rid != 0)
				rid_offsets[oreader.index.ent_meta[i].rid].range = rid_offsets[oreader.index.ent_meta[i].rid - 1].range;
		}
		rid_offsets[oreader.index.ent_meta[i].rid].min = oreader.index.ent_meta[i].minpos;
		rid_offsets[oreader.index.ent_meta[i].rid].max = oreader.index.ent_meta[i].maxpos;
	}
	std::cerr << "range=" << range << std::endl;
	for(int i = 0; i < rid_offsets.size(); ++i){
		std::cerr << "rid=" << i << "=" << rid_offsets[i].range << " -> " << rid_offsets[i].min << "-" << rid_offsets[i].max << std::endl;
	}

	// Step 2: Prepare n-tensor for storing output data.
	//         Matrix dimensions (1,2) correspond to pixel equivalents.
	//         Tensor dimensions (3,..) correspond to summary statistics for
	//         each bin (pixel).
	uint32_t xrange = std::ceil((float)range / x_bins);
	uint32_t yrange = std::ceil((float)range / y_bins);

	std::vector< std::vector<sstats> > mat(x_bins, std::vector<sstats>(y_bins));

	std::cerr << "partition range=" << range / x_bins << " and " << range / y_bins << std::endl;
	while(oreader.NextRecord()){
		if(settings.filter.Filter(oreader.it.rcd)){
			//writer.Add(*oreader.it.rcd);
			//std::cerr << oreader.it.rcd->Apos << "->" << oreader.it.rcd->Bpos << "\t" << (oreader.it.rcd->Apos/xrange) << "," << (oreader.it.rcd->Bpos/yrange) << "/" << x_bins << std::endl;
			if((oreader.it.rcd->Apos - rid_offsets[oreader.it.rcd->ridA].min)/xrange >= x_bins){
				std::cerr << "a=" << (oreader.it.rcd->Apos - rid_offsets[oreader.it.rcd->ridA].min) << "->" << (oreader.it.rcd->Apos - rid_offsets[oreader.it.rcd->ridA].min)/xrange << "/" << xrange << std::endl;
				std::cerr << oreader.it.rcd->Apos << "->" << oreader.it.rcd->Bpos << "\t" << (oreader.it.rcd->Apos/xrange) << "," << (oreader.it.rcd->Bpos/yrange) << "/" << x_bins << std::endl;
				exit(1);
			}
			if((oreader.it.rcd->Bpos - rid_offsets[oreader.it.rcd->ridB].min)/yrange >= y_bins){
				std::cerr << "b=" << (oreader.it.rcd->Bpos - rid_offsets[oreader.it.rcd->ridB].min) << "->" << (oreader.it.rcd->Bpos - rid_offsets[oreader.it.rcd->ridB].min)/range << "/" << yrange << std::endl;
				std::cerr << oreader.it.rcd->Apos << "->" << oreader.it.rcd->Bpos << "\t" << (oreader.it.rcd->Apos/xrange) << "," << (oreader.it.rcd->Bpos/yrange) << "/" << x_bins << std::endl;
				exit(1);
			}

			// Invoke summary statistics function.
			(mat[(oreader.it.rcd->Apos - rid_offsets[oreader.it.rcd->ridA].min)/xrange][(oreader.it.rcd->Bpos - rid_offsets[oreader.it.rcd->ridB].min)/yrange].*f)(oreader.it.rcd);
		}
	}
	std::cerr << "done" << std::endl;

	for(int i = 0; i < x_bins; ++i){
		std::cout << (mat[i][0].*r)(min_cutoff);
		for(int j = 1; j < y_bins; ++j){
			 std::cout << "\t" << (mat[i][j].*r)(min_cutoff);
		}
		std::cout << '\n';
	}
	std::cout.flush();
	std::cerr << "done" << std::endl;

	return 1;


	// todo: if data is sorted then only visit over*/lapping blocks
	if(settings.ivals.size()){
		std::cerr << settings.intervals.overlap_blocks.size() << std::endl;
		for(int i = 0; i < settings.intervals.overlap_blocks.size(); ++i){
			oreader.stream->seekg(settings.intervals.overlap_blocks[i]->foff);
			if(oreader.NextBlock() == false){
				std::cerr << "failed to get next block" << std::endl;
				return 1;
			}

			for(int j = 0; j < oreader.it.blk.n; ++j){
				assert(oreader.NextRecord());
				//oreader.it.rcd->PrintLD(std::cerr);
				if(settings.intervals.FilterInterval(*oreader.it.rcd)){
					continue;
				}

				if(settings.filter.Filter(oreader.it.rcd)){
					//writer.Add(*oreader.it.rcd);
				}
			}
		}
	} else if(settings.ivals.size()){
		while(oreader.NextRecord()){
			if(settings.intervals.FilterInterval(*oreader.it.rcd)){
				continue;
			}

			if(settings.filter.Filter(oreader.it.rcd)){
				//writer.Add(*oreader.it.rcd);
			}
		}
	} else {
		while(oreader.NextRecord()){
			if(settings.filter.Filter(oreader.it.rcd)){
				//writer.Add(*oreader.it.rcd);
			}
		}
	}

	return(0);
}
