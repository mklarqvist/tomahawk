#include <queue>

#include "two_reader.h"
#include "intervals.h"
#include "writer.h"
#include "two_sorter_structs.h"
#include "aggregation.h"

namespace tomahawk {

bool twk1_two_iterator::NextBlockRaw(){
	if(stream->good() == false){
		std::cerr << "stream died" << std::endl;
		return false;
	}

	//std::cerr << "pos=" << stream->tellg() << std::endl;

	uint8_t marker = 0;
	DeserializePrimitive(marker, *stream);
	if(marker == 0){
		//std::cerr << "0 marker found. stopping" << std::endl;
		return false;
	}
	//std::cerr << "marker=" << (int)marker << std::endl;
	//assert(marker == 1);
	if(marker != 1){
		std::cerr << "Marker!=1 is " << (int)marker << " @ " << stream->tellg() << " good=" << stream->good() << std::endl;
		exit(1);
	}

	*stream >> oblk;
	if(stream->good() == false){
		std::cerr << "stream died" << std::endl;
		return false;
	}

	assert(oblk.bytes.size() == oblk.nc);
	buf.resize(oblk.n);
	offset = 0;
	rcd = nullptr;

	return true;
}

bool twk1_two_iterator::NextBlock(){
	if(this->NextBlockRaw() == false)
		return false;

	// Decompress data
	zcodec.Decompress(oblk.bytes, buf);
	buf >> blk;
	buf.reset();
	if(blk.n) rcd = &blk.rcds[0];

	return true;
}

bool twk1_two_iterator::NextRecord(){
	if(offset == blk.n){
		if(this->NextBlock() == false)
			return false;

		offset = 0;
	}
	rcd = &blk.rcds[offset++];
	return true;
}

two_reader::two_reader() : buf(nullptr), stream(nullptr){}
two_reader::~two_reader(){ delete stream; }

bool two_reader::BuildIntervals(std::vector<std::string>& strings, const uint32_t n_contigs,
		           const IndexOutput& index, const VcfHeader& hdr)
{
	return(intervals.Build(strings, n_contigs, index, hdr));
}

bool two_reader::FilterInterval(const twk1_two_t* rec) const { return(intervals.FilterInterval(*rec)); }
bool two_reader::FilterInterval(const twk1_two_t& rec) const { return(intervals.FilterInterval(rec)); }

IndexEntryOutput* two_reader::GetIntervalBlock(const uint32_t p){ return(intervals.GetOverlapBlock(p)); }
const std::vector<IndexEntryOutput*>& two_reader::GetIntervalBlocks() const { return(intervals.overlap_blocks); }

bool two_reader::Open(std::string file){
	fstream.open(file, std::ios::in|std::ios::binary|std::ios::ate);
	if(!fstream.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to open \"" << file << "\"..." << std::endl;
		return false;
	}
	buf = fstream.rdbuf();
	stream = new std::istream(buf);

	uint64_t filesize = stream->tellg();
	stream->seekg(0);

	// read magic
	char magic[TOMAHAWK_LD_MAGIC_HEADER_LENGTH];
	stream->read(magic, TOMAHAWK_LD_MAGIC_HEADER_LENGTH);
	if(strncmp(magic, TOMAHAWK_LD_MAGIC_HEADER.data(), TOMAHAWK_LD_MAGIC_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR") << "Failed to read TWO magic string! Incorrect file format!" << std::endl;
		return false;
	}

	// Read, decompress, and parse header
	uint64_t buf_size = 0, obuf_size = 0;
	stream->read(reinterpret_cast<char*>(&buf_size), sizeof(uint64_t));
	stream->read(reinterpret_cast<char*>(&obuf_size),sizeof(uint64_t));
	twk_buffer_t obuf(obuf_size);
	twk_buffer_t buf(buf_size);
	stream->read(obuf.data(),obuf_size);
	obuf.n_chars_ = obuf_size;

	if(zcodec.Decompress(obuf, buf) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to decompress header! Corrupted file!" << std::endl;
		return false;
	}
	assert(buf.size() == buf_size);
	buf >> hdr;
	buf.reset(); obuf.reset();

	// Remember seek point to start of data.
	uint64_t data_start = stream->tellg();

	// seek to end-of-file
	// seek back to end of file marker and position where index offset is stored
	stream->seekg(filesize - TOMAHAWK_FILE_EOF_LENGTH - sizeof(uint64_t));
	uint64_t offset_start_index = 0;
	stream->read(reinterpret_cast<char*>(&offset_start_index), sizeof(uint64_t));

	// Seek to start of offst
	stream->seekg(offset_start_index);
	if(stream->good() == false){
		std::cerr << utility::timestamp("ERROR") << "Failed seek in file! Corrupted file!" << std::endl;
		return false;
	}

	// Load index
	uint8_t marker = 0;
	stream->read(reinterpret_cast<char*>(&marker),   sizeof(uint8_t));
	stream->read(reinterpret_cast<char*>(&buf_size), sizeof(uint64_t));
	stream->read(reinterpret_cast<char*>(&obuf_size),sizeof(uint64_t));
	obuf.resize(obuf_size), buf.resize(buf_size);
	stream->read(obuf.data(),obuf_size);
	obuf.n_chars_ = obuf_size;

	if(zcodec.Decompress(obuf, buf) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to decompress! Corrupted file!" << std::endl;
		return false;
	}
	buf >> index;

	// Seek back to the beginning of data.
	stream->seekg(data_start);

	// Assign stream.
	it.stream = stream;

	return(stream->good());
}

bool two_reader::Sort(){
	two_sorter_settings settings;
	//this->settings = settings;
	return(Sort(settings));
}

bool two_reader::Sort(two_sorter_settings& settings){
	if(settings.in.length() == 0){
		std::cerr << utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return false;
	}

	// File reader.
	if(Open(settings.in) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to open \"" << settings.in << "\"..."  << std::endl;
		return false;
	}

	twk1_two_block_t blk2;
	twk_buffer_t obuf, obuf2;

	// Distrubution.
	uint64_t b_unc = 0, n_recs = 0;
	std::cerr << utility::timestamp("LOG") << "Blocks: " << utility::ToPrettyString(index.n) << std::endl;
	for(int i = 0; i < index.n; ++i){
		b_unc  += index.ent[i].b_unc;
		n_recs += index.ent[i].n;
	}
	std::cerr << utility::timestamp("LOG") << "Uncompressed size: " << utility::ToPrettyDiskString(b_unc) << std::endl;
	std::cerr << utility::timestamp("LOG") << "Sorting " << utility::ToPrettyString(n_recs) << " records..." << std::endl;

	if(b_unc == 0){
		std::cerr << utility::timestamp("ERROR") << "Cannot sort empty file..." << std::endl;
		return false;
	}

	if(index.n < settings.n_threads) settings.n_threads = index.n;
	uint64_t b_unc_thread = b_unc / settings.n_threads;
	std::cerr << utility::timestamp("LOG","THREAD") << "Data/thread: " << utility::ToPrettyDiskString(b_unc_thread) << std::endl;

	std::vector< std::pair<uint32_t,uint32_t> > ranges;
	uint64_t f = 0, t = 0, b_unc_tot = 0;
	for(int i = 0; i < index.n; ++i){
		if(b_unc_tot >= b_unc_thread){
			ranges.push_back(std::pair<uint32_t,uint32_t>(f, t));
			b_unc_tot = 0;
			f = t;
		}
		b_unc_tot += index.ent[i].b_unc;
		//std::cerr << "index=" << index.ent[i].b_cmp << " -> " << index.ent[i].b_unc << std::endl;
		++t;
	}
	if(f != t){
		ranges.push_back(std::pair<uint32_t,uint32_t>(f, t));
		b_unc_tot = 0;
		f = t;
	}
	//std::cerr << "ranges: " << ranges.size() << std::endl;
	assert(ranges.back().second == index.n);
	assert(ranges.size() <= settings.n_threads);

	//return(1);

	twk_sort_progress progress_sort;
	progress_sort.n_cmps = n_recs;
	std::thread* psthread = progress_sort.Start();

	twk_sort_slave* slaves = new twk_sort_slave[settings.n_threads];
	uint32_t range_thread = index.n / settings.n_threads;
	for(int i = 0; i < settings.n_threads; ++i){
		slaves[i].f = ranges[i].first;
		slaves[i].t = ranges[i].second;
		slaves[i].m_limit = settings.memory_limit;
		slaves[i].filename = settings.in;
		slaves[i].c_level = settings.c_level;
		slaves[i].progress = &progress_sort;

		std::string suffix    = twk_two_writer_t::RandomSuffix();
		std::string base_path = twk_two_writer_t::GetBasePath(settings.out);
		std::string base_name = twk_two_writer_t::GetBaseName(settings.out);
		std::string temp_out  = (base_path.size() ? base_path + "/" : "") + base_name + "_" + suffix + ".two";
		slaves[i].tmp_filename = temp_out;
		std::cerr << utility::timestamp("LOG","THREAD") << "Slave-" << i << ": range=" << slaves[i].f << "->" << slaves[i].t << "/" << index.n << " and name " << slaves[i].tmp_filename << std::endl;
	}

	for(int i = 0; i < settings.n_threads; ++i){
		if(slaves[i].Start(index) == nullptr){
			std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to spawn slave" << std::endl;
			return false;
		}
	}
	for(int i = 0; i < settings.n_threads; ++i) slaves[i].thread->join();
	progress_sort.is_ticking = false;
	progress_sort.PrintFinal();

	// temp
	/*
	for(int i = 0; i < settings.n_threads; ++i){
		std::cerr << i << "\t" << slaves[i].run_ivals.size() << std::endl;
		for(int j = 0; j < slaves[i].run_ivals.size(); ++j){
			std::cerr << "\tblock-" << j << ": " << slaves[i].run_ivals[j].size() << "\t(" << slaves[i].run_ivals[j][0].ref_rid << "," << slaves[i].run_ivals[j][0].n_run << "," << slaves[i].run_ivals[j][0].minp << "-" << slaves[i].run_ivals[j][0].maxp << ")";
			for(int k = 1; k < slaves[i].run_ivals[j].size(); ++k){
				std::cerr << ", (" << slaves[i].run_ivals[j][k].ref_rid << "," << slaves[i].run_ivals[j][k].n_run << "," << slaves[i].run_ivals[j][k].minp << "-" << slaves[i].run_ivals[j][k].maxp << ")";
			}
			std::cerr << std::endl;
		}
	}
	*/
	//

	uint32_t n_queues = 0;
	for(int i = 0; i < settings.n_threads; ++i){
		for(int j = 0; j < slaves[i].local_idx.size(); ++j){
			++n_queues;
		}
	}


	// Merge
	obuf.reset(); obuf2.reset();
	obuf.resize(256000);
	obuf2.resize(256000);

	//uint32_t k = 0;
	uint64_t maxmem_queue = settings.memory_limit * settings.n_threads * 1e9 / 15; // assume compression ratio is 15
	uint64_t mem_queue = maxmem_queue / n_queues;
	mem_queue = mem_queue < sizeof(twk1_two_t) ? sizeof(twk1_two_t) : mem_queue;

	std::priority_queue<two_queue_entry> queue;

	std::cerr << utility::timestamp("LOG") << "Spawning " << utility::ToPrettyString(n_queues) << " queues with " << utility::ToPrettyDiskString(mem_queue) << " each..." << std::endl;
	twk_two_stream_iterator* its = new twk_two_stream_iterator[n_queues];
	twk1_two_t rec;
	uint64_t n_rec_total = 0;
	uint32_t local_queue = 0;
	for(int i = 0; i < settings.n_threads; ++i){
		for(int j = 0; j < slaves[i].local_idx.size(); ++j, ++local_queue){
			// open iterators
			if(its[local_queue].Open(slaves[i].tmp_filename,
									 slaves[i].local_idx[j].foff,
									 slaves[i].local_idx[j].fend,
									 slaves[i].local_idx[j].n,
									 slaves[i].local_idx[j].nc) == false)
			{
				std::cerr << utility::timestamp("ERROR") << "Failed open \"" << slaves[i].tmp_filename << "\"..." << std::endl;
				return false;
			}

			if(its[local_queue].Next(rec, mem_queue) == false){
				std::cerr << utility::timestamp("ERROR") << "Failed to get next" << std::endl;
				return false;
			}

			queue.push(two_queue_entry(rec, local_queue));

			n_rec_total += slaves[i].local_idx[j].n / twk1_two_t::packed_size;
		}
	}

	if(queue.empty()){
		std::cerr << utility::timestamp("ERROR","SORT") << "No data in queue..." << std::endl;
		return false;
	}

	twk_two_writer_t owriter;
	owriter.oindex.SetChroms(hdr.GetNumberContigs());
	if(settings.out.size() == 0 || (settings.out.size() == 1 && settings.out == "-")){
		std::cerr << utility::timestamp("LOG","WRITER") << "Writing to stdout..." << std::endl;
	} else {
		std::string extension = twk_two_writer_t::GetExtension(settings.out);
		if(extension != "two"){
			settings.out += ".two";
		}
		std::cerr << utility::timestamp("LOG","WRITER") << "Opening \"" << settings.out << "\"..." << std::endl;
	}

	if(owriter.Open(settings.out) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed top open \"" << settings.out << "\"..." << std::endl;
		return false;
	}
	owriter.mode = 'b';
	owriter.oindex.state = TWK_IDX_SORTED;
	owriter.SetCompressionLevel(settings.c_level);
	// Write header
	std::string sort_string = "\n##tomahawk_sortVersion=" + std::string(VERSION) + "\n";
	sort_string += "##tomahawk_sortCommand=" + LITERAL_COMMAND_LINE + "; Date=" + utility::datetime() + "\n";
	hdr.literals_ += sort_string;
	if(owriter.WriteHeader(*this) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to write header!" << std::endl;
		return false;
	}

	// Reference
	uint32_t ridA = queue.top().rec.ridA;

	Timer timer; timer.Start();
	uint64_t n_entries_out = 0;
	twk_sort_progress progress;
	progress.n_cmps = n_recs;
	std::thread* pthread = progress.Start();

	while(queue.empty() == false){
		// peek at top entry in queue
		const uint32_t id = queue.top().qid;

		if(queue.top().rec.ridA != ridA){
			if(owriter.WriteBlock() == false){
				std::cerr << utility::timestamp("ERROR") << "Failed to flush block..." << std::endl;
				return false;
			}
		}
		owriter.Add(queue.top().rec);
		ridA = queue.top().rec.ridA;

		++progress.cmps;

		// remove this record from the queue
		queue.pop();

		while(its[id].Next(rec, mem_queue)){
			if(!(rec < queue.top().rec)){
				queue.push( two_queue_entry(rec, id) );
				break;
			}

			if(rec.ridA != ridA){
				if(owriter.WriteBlock() == false){
					std::cerr << utility::timestamp("ERROR") << "Failed to flush block..." << std::endl;
					return false;
				}
			}
			owriter.Add(rec);
			ridA = rec.ridA;

			++progress.cmps;
		}
	}
	progress.is_ticking = false;
	progress.PrintFinal();

	owriter.flush();
	owriter.WriteFinal();
	owriter.close();
	std::cerr << utility::timestamp("LOG") << "Finished merging! Time: " << timer.ElapsedString() << std::endl;
	//std::cerr << "deleting intermediary" << std::endl;

	std::cerr << utility::timestamp("LOG") << "Deleting temp files..." << std::endl;
	std::cerr.flush();
	for(int i = 0; i < settings.n_threads; ++i){
	if( remove( slaves[i].tmp_filename.c_str() ) != 0 ){
		std::cerr << utility::timestamp("ERROR") << "Error deleting file " << slaves[i].tmp_filename << "!" << std::endl;
	} else {
		std::cerr << utility::timestamp("LOG") << "Deleted " << slaves[i].tmp_filename << std::endl;
	  }
	}

	delete[] slaves;
	delete[] its;
	std::cerr << utility::timestamp("LOG") << "Finished!" << std::endl;
	return true;
}

bool two_reader::Decay(twk_two_settings& settings, int64_t window_bp, int32_t n_bins){
	if(window_bp <= 0){
		std::cerr << utility::timestamp("ERROR") << "Window size cannot be <= 0 (provided " << window_bp << ")..." << std::endl;
		return false;
	}

	if(n_bins <= 0){
		std::cerr << utility::timestamp("ERROR") << "Number of bins cannot be <= 0 (provided " << n_bins << ")..." << std::endl;
		return false;
	}

	if(settings.in.length() == 0){
		std::cerr << utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return false;
	}

	// File reader.
	if(Open(settings.in) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to open \"" << settings.in << "\"..."  << std::endl;
		return false;
	}

	// Build intervals data structures if any are available.
	if(BuildIntervals(settings.ivals, hdr.GetNumberContigs(),
	                   index, hdr) == false)
	{
		return false;
	}

	const uint32_t n_bins_i = n_bins - 1;
	uint32_t n_range_bin = window_bp / n_bins;
	std::vector< std::pair<double, uint64_t> > decay(n_bins, {0,0});

	while(NextRecord()){
		// Same contig only.
		if(it.rcd->ridA == it.rcd->ridB){
			// Upper trig only.
			if(it.rcd->Apos < it.rcd->Bpos){
				decay[std::min((it.rcd->Bpos - it.rcd->Apos) / n_range_bin, n_bins_i)].first += it.rcd->R2;
				++decay[std::min((it.rcd->Bpos - it.rcd->Apos) / n_range_bin, n_bins_i)].second;
			}
		}
	}

	std::cout << "From\tTo\tMean\tFrequency\n";
	for(int i = 0; i < decay.size(); ++i){
		std::cout << (i*n_range_bin) << '\t' << ((i+1)*n_range_bin) << '\t' << decay[i].first/std::max(decay[i].second,(uint64_t)1) << '\t' << decay[i].second << '\n';
	}
	std::cout.flush();

	return true;
}

bool two_reader::PositionalDecay(twk_two_settings& settings){
    if(settings.in.length() == 0){
        std::cerr << utility::timestamp("ERROR") << "No input value specified..." << std::endl;
        return false;
    }

    // File reader.
    if(Open(settings.in) == false){
        std::cerr << utility::timestamp("ERROR") << "Failed to open \"" << settings.in << "\"..."  << std::endl;
        return false;
    }

    // Build intervals data structures if any are available.
    if(BuildIntervals(settings.ivals, hdr.GetNumberContigs(),
                       index, hdr) == false)
    {
        return false;
    }

    if(NextRecord() == false){
        std::cerr << utility::timestamp("ERROR") << "Failed to get a record..."  << std::endl;
        return false;
    }

    //
    struct twk_sstats_pos : public twk_sstats {
        twk_sstats_pos() : rid(0), pos(0){}
        twk_sstats_pos(uint32_t chrom, uint32_t position) : rid(chrom), pos(position){}

        void operator+=(const twk1_two_t* rec){
            if(rec->ridA == rec->ridB && rec->Apos < rec->Bpos)
                Add(rec->Bpos, 1);
        }

        uint32_t rid, pos;
    };
    //

    uint32_t rid_prev = it.rcd->ridA;
    uint32_t pos_prev = it.rcd->Apos;
    std::vector<twk_sstats_pos> variants;
    variants.push_back(twk_sstats_pos(it.rcd->ridA, it.rcd->Apos));

    while(NextRecord()){
        if(FilterInterval(it.rcd)) continue;

        if(it.rcd->ridA != rid_prev || it.rcd->Apos != pos_prev){
            variants.push_back(twk_sstats_pos(it.rcd->ridA, it.rcd->Apos));
            rid_prev = it.rcd->ridA;
            pos_prev = it.rcd->Apos;
            variants.back().n = 1;
        }

        variants.back() += it.rcd;
    }

    std::cerr << "variants = " << variants.size() << std::endl;
    std::cout << std::fixed;
    for(int i = 0; i < variants.size(); ++i){
        std::cout << variants[i].rid << '\t' << variants[i].pos << '\t' << variants[i].n << '\t' << variants[i].GetMean()-variants[i].pos << '\n';
    }
    std::cout.flush();

    return true;
}

// aggregation
bool two_reader::Aggregate(twk1_aggregate_t& agg,
		twk_two_settings& settings,
		std::string agg_name, std::string red_name,
		int32_t xbins, int32_t ybins,
		int32_t min_cutoff,
		bool verbose, bool progress)
{
	if(agg_name.size() == 0){
		if(verbose) std::cerr << utility::timestamp("ERROR") << "No aggregation function (-f) provided..." << std::endl;
		return(false);
	}

	if(red_name.size() == 0){
		if(verbose) std::cerr << utility::timestamp("ERROR") << "No reduce function (-r) provided..." << std::endl;
		return(false);
	}

	if(min_cutoff < 0){
		if(verbose) std::cerr << utility::timestamp("ERROR") << "Cannot have a min-cutoff (-c) < 0..." << std::endl;
		return(false);
	}

	if(settings.n_threads <= 0){
		std::cerr << utility::timestamp("ERROR") << "Cannot have <= 0 threads (-t)..." << std::endl;
		return(false);
	}
	twk_sstats::aggfunc f = &twk_sstats::AddR2;
	twk_sstats::redfunc r = &twk_sstats::GetMean;

	// Transform string of aggregation function name into lower then try to map
	// name to existing function names.
	std::transform(agg_name.begin(), agg_name.end(), agg_name.begin(), ::tolower);
	if(agg_name == "r2")         { f = &twk_sstats::AddR2;   }
	else if(agg_name == "r")     { f = &twk_sstats::AddR;    }
	else if(agg_name == "d")     { f = &twk_sstats::AddD;    }
	else if(agg_name == "dprime"){ f = &twk_sstats::AddDprime; }
	else if(agg_name == "dp")    { f = &twk_sstats::AddDprime; }
	else if(agg_name == "p")     { f = &twk_sstats::AddP;    }
	else if(agg_name == "hets")  { f = &twk_sstats::AddHets; }
	else if(agg_name == "alts")  { f = &twk_sstats::AddAlts; }
	else if(agg_name == "het")   { f = &twk_sstats::AddHets; }
	else if(agg_name == "alt")   { f = &twk_sstats::AddAlts; }
	else {
		std::cerr << utility::timestamp("ERROR") << "Unknown aggregation function \"" << agg_name << "\"..." << std::endl;
		return(false);
	}

	// Transform string of reduction function name into lower then try to map
	// name to existing function names.
	std::transform(red_name.begin(), red_name.end(), red_name.begin(), ::tolower);
	if(red_name == "mean")          { r = &twk_sstats::GetMean;  }
		else if(red_name == "max")  { r = &twk_sstats::GetMax;   }
		else if(red_name == "min")  { r = &twk_sstats::GetMin;   }
		else if(red_name == "count"){ r = &twk_sstats::GetCount; }
		else if(red_name == "n")    { r = &twk_sstats::GetCount; }
		else if(red_name == "total"){ r = &twk_sstats::GetTotal; }
		else if(red_name == "sd")   { r = &twk_sstats::GetStandardDeviation; }
		else {
			std::cerr << utility::timestamp("ERROR") << "Unknown reduce function \"" << red_name << "\"..." << std::endl;
			return(false);
		}

	if(xbins < 5){
		std::cerr << utility::timestamp("ERROR") << "Number of x-bins cannot be < 5!" << std::endl;
		return(false);
	}

	if(ybins < 5){
		std::cerr << utility::timestamp("ERROR") << "Number of y-bins cannot be < 5!" << std::endl;
		return(false);
	}

	if(settings.in.length() == 0){
		std::cerr << utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(false);
	}

	// Open file handle.
	if(Open(settings.in) == false) return 1;


	// Build intervals data structures if any are available.
	/*if(settings.intervals.Build(settings.ivals,
								hdr.GetNumberContigs(),
								index,
								hdr) == false)
	{
		return 1;
	}*/

	// Construct filters.
	//settings.filter.Build();

	// Print messages
	if(verbose){
		ProgramMessage();
		std::cerr << utility::timestamp("LOG") << "Calling aggregate..." << std::endl;
	}

	// Algorithmic overview.
	// Step 1: Find maximum and minimum X and Y values.
	// Step 2: Partition into (maxY-minY)/#bins and (maxX-minX)/#bins buckets.
	//      a: Range can be either dynamic min and max given the data, or;
	//      b: Interval (from,to)-tuple.
	// Step 3: Iterate over data and update summary statistics in buckets.
	// Step 4: Output data.
	if(verbose) std::cerr << utility::timestamp("LOG") << "Performing 2-pass over data..." << std::endl;

	// Step 1: iterate over index entries and find what contigs are used
	if(verbose) std::cerr << utility::timestamp("LOG") << "===== First pass (peeking at landscape) =====" << std::endl;
	// Distrubution.
	uint64_t b_unc = 0, n_recs = 0;
	if(verbose) std::cerr << utility::timestamp("LOG") << "Blocks: " << utility::ToPrettyString(index.n) << std::endl;
	for(int i = 0; i < index.n; ++i){
		b_unc  += index.ent[i].b_unc;
		n_recs += index.ent[i].n;
	}
	if(verbose) std::cerr << utility::timestamp("LOG") << "Uncompressed size: " << utility::ToPrettyDiskString(b_unc) << std::endl;

	if(b_unc == 0){
		if(verbose) std::cerr << utility::timestamp("LOG") << "Cannot aggregate empty file..." << std::endl;
		return false;
	}

	if(index.n < settings.n_threads) settings.n_threads = index.n;
	uint64_t b_unc_thread = b_unc / settings.n_threads;
	if(verbose) std::cerr << utility::timestamp("LOG","THREAD") << "Data/thread: " << utility::ToPrettyDiskString(b_unc_thread) << std::endl;

	std::vector< std::pair<uint32_t,uint32_t> > ranges;
	uint64_t fR = 0, tR = 0, b_unc_tot = 0;
	for(int i = 0; i < index.n; ++i){
		if(b_unc_tot >= b_unc_thread){
			ranges.push_back(std::pair<uint32_t,uint32_t>(fR, tR));
			b_unc_tot = 0;
			fR = tR;
		}
		b_unc_tot += index.ent[i].b_unc;
		++tR;
	}
	if(fR != tR){
		ranges.push_back(std::pair<uint32_t,uint32_t>(fR, tR));
		b_unc_tot = 0;
		fR = tR;
	}
	assert(ranges.back().second == index.n);
	assert(ranges.size() <= settings.n_threads);

	twk_sort_progress progress_sort;
	progress_sort.n_cmps = n_recs;
	if(progress) progress_sort.Start();

	twk_agg_slave* slaves = new twk_agg_slave[settings.n_threads];
	uint32_t range_thread = index.n / settings.n_threads;
	for(int i = 0; i < settings.n_threads; ++i){
		slaves[i].f = ranges[i].first;
		slaves[i].t = ranges[i].second;
		slaves[i].filename = settings.in;
		slaves[i].progress = &progress_sort;
		//std::cerr << "thread-" << i << " " << slaves[i].f << "-" << slaves[i].t << std::endl;
	}

	for(int i = 0; i < settings.n_threads; ++i){
		if(slaves[i].StartFindRanges(*this) == nullptr){
			if(verbose) std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to spawn slave" << std::endl;
			return false;
		}
	}
	for(int i = 0; i < settings.n_threads; ++i) slaves[i].thread->join();
	progress_sort.is_ticking = false;
	if(progress && verbose) progress_sort.PrintFinal();

	// Reduce.
	for(int i = 1; i < settings.n_threads; ++i){
		for(int j = 0; j < slaves[0].contig_avail.size(); ++j){
			slaves[0].contig_avail[j].set = std::max(slaves[0].contig_avail[j].set, slaves[i].contig_avail[j].set);
			slaves[0].contig_avail[j].min = std::min(slaves[0].contig_avail[j].min, slaves[i].contig_avail[j].min);
			slaves[0].contig_avail[j].max = std::max(slaves[0].contig_avail[j].max, slaves[i].contig_avail[j].max);
		}
	}

	// Reduce.
	uint32_t n_chrom_set = slaves[0].contig_avail[0].set;
	for(int i = 0; i < slaves[0].contig_avail.size(); ++i){
		n_chrom_set += slaves[0].contig_avail[i].set;
	}
	//std::cerr << "chrom set=" << n_chrom_set << std::endl;

	// Step 2: Determine boundaries given the contigs that were set.
	//         Calculate the landscape ranges (X and Y dimensions).
	// Approach 2: Dropping regions with no data.
	uint64_t range = 0;
	std::vector<twk1_aggregate_t::offset_tuple> rid_offsets(slaves[0].contig_avail.size());

	/**<
	 * If there is only chromosome set then restrict the (X,Y) landscape to the
	 * available data range. This is in contrast to cases where N > 1, where we
	 * consider the entire genomic range of the affected chromosomes irrespective
	 * of how much range is actually used.
	 */
	if(n_chrom_set == 1){
		if(slaves[0].contig_avail[0].set){
			rid_offsets[0].range = slaves[0].contig_avail[0].max - slaves[0].contig_avail[0].min + 1;
			range += slaves[0].contig_avail[0].max - slaves[0].contig_avail[0].min + 1;
		} else
			rid_offsets[0].range = 0;

		rid_offsets[0].min = slaves[0].contig_avail[0].min;
		rid_offsets[0].max = slaves[0].contig_avail[0].max;

		for(int i = 1; i < slaves[0].contig_avail.size(); ++i){
			if(slaves[0].contig_avail[i].set){
				// Cumulative offset for the current rid equals the previous rid
				rid_offsets[i].range = rid_offsets[i - 1].range + (slaves[0].contig_avail[i].max - slaves[0].contig_avail[i].min + 1);
				range += (slaves[0].contig_avail[i].max - slaves[0].contig_avail[i].min + 1);
			} else {
				rid_offsets[i].range = rid_offsets[i - 1].range;
			}
			rid_offsets[i].min = slaves[0].contig_avail[i].min;
			rid_offsets[i].max = slaves[0].contig_avail[i].max;
		}

	}
	// If there is data from n>1 chromosomes.
	else {
		if(slaves[0].contig_avail[0].set){
			rid_offsets[0].range = hdr.contigs_[0].n_bases;
			range += hdr.contigs_[0].n_bases;
		} else
			rid_offsets[0].range = 0;

		rid_offsets[0].min = 0;
		rid_offsets[0].max = hdr.contigs_[0].n_bases;

		for(int i = 1; i < slaves[0].contig_avail.size(); ++i){
			if(slaves[0].contig_avail[i].set){
				// Cumulative offset for the current rid equals the previous rid
				rid_offsets[i].range = rid_offsets[i - 1].range + hdr.contigs_[i].n_bases;
				range += hdr.contigs_[i].n_bases;

			} else {
				rid_offsets[i].range = rid_offsets[i - 1].range;
			}
			rid_offsets[i].min = 0;
			rid_offsets[i].max = hdr.contigs_[i].n_bases;
		}
	}

	//std::cerr << "range=" << range << std::endl;
	//for(int i = 0; i < rid_offsets.size(); ++i){
	//	std::cerr << "rid=" << i << "=" << rid_offsets[i].range << " -> " << rid_offsets[i].min << "-" << rid_offsets[i].max << std::endl;
	//}

	// Step 3: Second pass over data.
	//         Prepare n-tensor for storing output data.
	//         Matrix dimensions (1,2) correspond to pixel equivalents.
	//         Tensor dimensions (3,..) correspond to summary statistics for
	//         each bin (pixel).
	uint32_t xrange = std::ceil((float)range / xbins);
	uint32_t yrange = std::ceil((float)range / ybins);

	//twk1_aggregate_t agg(xbins, ybins);
	agg.n = xbins*ybins;
	agg.x = xbins;
	agg.y = ybins;
	agg.data = new double[agg.n];
	agg.bpx = xrange;
	agg.bpy = yrange;
	agg.range = range;
	agg.rid_offsets = rid_offsets;
	agg.n_original = n_recs;

	if(verbose){
		std::cerr << utility::timestamp("LOG") << "===== Second pass (building matrix) =====" << std::endl;
		std::cerr << utility::timestamp("LOG") << "Aggregating " << utility::ToPrettyString(n_recs) << " records..." << std::endl;
		std::cerr << utility::timestamp("LOG","THREAD") << "Allocating: " << utility::ToPrettyDiskString(sizeof(twk_sstats)*xbins*ybins*settings.n_threads) << " for matrices..." << std::endl;
	}

	twk_sort_progress progress_sort_step2;
	progress_sort_step2.n_cmps = n_recs;
	if(progress) progress_sort_step2.Start();
	//std::cerr << "range=" << range << " x,y = " << xrange << " bp/pixel " << " and " << yrange << " bp/pixel" << std::endl;
	for(int i = 0; i < settings.n_threads; ++i){
		slaves[i].rid_offsets = rid_offsets;
		slaves[i].progress = &progress_sort_step2;
		if(slaves[i].StartBuildMatrix(*this, f, r, xbins, ybins, xrange, yrange) == nullptr){
			if(verbose) std::cerr << utility::timestamp("ERROR","THREAD") << "Failed to spawn slave" << std::endl;
			return false;
		}
	}
	for(int i = 0; i < settings.n_threads; ++i) slaves[i].thread->join();
	progress_sort_step2.is_ticking = false;
	if(progress && verbose) progress_sort_step2.PrintFinal();
	for(int i = 1; i < settings.n_threads; ++i) slaves[0].AddMatrix(slaves[i]);

	// Print matrix
	//slaves[0].PrintMatrix(std::cout, min_cutoff);
	slaves[0].Overload(agg, min_cutoff);
	//std::cout << agg;

	if(verbose){
		std::cerr << utility::timestamp("LOG") << "Aggregated " << utility::ToPrettyString(n_recs) << " records in " << utility::ToPrettyString(xbins*ybins) << " bins." << std::endl;
		std::cerr << utility::timestamp("LOG") << "Finished." << std::endl;
	}

	delete[] slaves;
	return(true);
}


}
