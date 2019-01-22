#include "intervals.h"

namespace tomahawk {

interval_pair_payload::interval_pair_payload() : rid(-1), offset(0), mate(0){}
interval_pair_payload::interval_pair_payload(const int32_t rid, const uint32_t mate_offset, const uint8_t type) :
	mate(type), rid(rid), offset(mate_offset)
{}

bool interval_pair_payload::operator<(const interval_pair_payload& other) const {
	if(rid < other.rid) return true;
	if(other.rid < rid) return false;
	if(offset < other.offset) return true;
	if(other.offset < offset) return false;
	return true;
}

twk_intervals::twk_intervals() : n_c(0), itree(nullptr){}
twk_intervals::twk_intervals(const uint32_t n_contigs) :
	n_c(n_contigs),
	itree(new algorithm::IntervalTree<uint32_t,uint32_t>*[n_c])
{
	// Reserve memory for 250 intervals for each rid.
	// This is to make sure that pointers of linked reads
	// doesn't become invalid when resizing these vectors.
	for(int i = 0; i < n_c; ++i) ivecs[i].reserve(250);
}

twk_intervals::~twk_intervals(){
	for(uint32_t i = 0; i < n_c; ++i) delete itree[i];
	delete[] itree;
}

void twk_intervals::Dedupe(){
	for(int i = 0; i < n_c; ++i){ // every rid
		if(ivecs[i].size() == 0) continue;
		std::vector<interval> ivals;
		std::sort(ivecs[i].begin(), ivecs[i].end());
		ivals.push_back(ivecs[i][0]);
		for(int j = 1; j < ivecs[i].size(); ++j){
			// Not-inclusive right-interval [from, to)
			if(  ivecs[i][j].start <  ivals.back().stop
			  && ivecs[i][j].stop  >= ivals.back().start)
			{
				ivals.back().stop = ivecs[i][j].stop;
			}
			else ivals.push_back(ivecs[i][j]);
		}
		ivecs[i] = ivals;
	}
}


bool twk_intervals::Build(const uint32_t n_contigs, const Index& index){
	for(uint32_t i = 0; i < n_c; ++i) delete itree[i];
	delete[] itree; itree = nullptr;

	if(n_contigs == 0) return false;
	n_c = n_contigs;
	this->Dedupe();
	for(uint32_t i = 0; i < ivecs.size(); ++i){
		for(int j = 0; j < ivecs[i].size(); ++j){
			std::vector<IndexEntry*> idx = index.FindOverlap(i,ivecs[i][j].start,ivecs[i][j].stop);
			overlap_blocks.insert(overlap_blocks.end(), idx.begin(), idx.end());
		}
	}

	itree = new algorithm::IntervalTree<uint32_t,uint32_t>*[n_c];
	for(uint32_t i = 0; i < n_c; ++i){
		itree[i] = new algorithm::IntervalTree<uint32_t,uint32_t>(ivecs[i]);
	}

	if(overlap_blocks.size() == 0){
		std::cerr << utility::timestamp("ERROR","INTERVAL") << "Found no blocks overlapping the provided range(s)..." << std::endl;
		return(false);
	}

	return(true);
}

bool twk_intervals::ParseIntervalStrings(const std::vector<std::string>& ivals, VcfHeader& hdr){
	if(ivals.size() == 0)
		return true;

	for(int i = 0; i < ivals.size(); ++i){
		if(this->ParseIntervalString(ivals[i], hdr) == false)
			return false;
	}
	return true;
}

bool twk_intervals::ParseIntervalString(const std::string& s, const VcfHeader& hdr){
	if(std::regex_match(s, TWK_REGEX_CONTIG_RANGE)){ // contig only
		std::vector<std::string> rid = utility::split(s, ':');
		std::vector<std::string> pos;
		if(rid.size() == 2){
			pos = utility::split(rid[1], '-');
		} else {
			std::cerr << utility::timestamp("ERROR","INTERVAL") << "Illegal format: " << s << std::endl;
			return false;
		}
		uint32_t ival_from = (uint32_t)std::atof(pos[0].c_str());
		uint32_t ival_to   = (uint32_t)std::atof(pos[1].c_str());
		const VcfContig* contig = hdr.GetContig(rid[0]);
		if(contig == nullptr){
			std::cerr << utility::timestamp("ERROR","INTERVAL") << "Contig does not exist in string " << s << std::endl;
			return false;
		}
		uint32_t ival_rid = contig->idx;
		ivecs[ival_rid].push_back(algorithm::Interval<uint32_t,uint32_t>(ival_from,ival_to,0));

		return true;

	} else if(std::regex_match(s, TWK_REGEX_CONTIG_POSITION)){ // contig and position
		std::vector<std::string> rid = utility::split(s, ':');

		uint32_t ival_from = (uint32_t)std::atof(rid[1].c_str());
		uint32_t ival_to   = (uint32_t)std::atof(rid[1].c_str()) + 1;
		const VcfContig* contig = hdr.GetContig(rid[0]);
		if(contig == nullptr){
			std::cerr << utility::timestamp("ERROR","INTERVAL") << "Contig does not exist in string " << s << std::endl;
			return false;
		}
		uint32_t ival_rid = contig->idx;
		ivecs[ival_rid].push_back(algorithm::Interval<uint32_t,uint32_t>(ival_from,ival_to,0));

		return true;
	} else if(std::regex_match(s, TWK_REGEX_CONTIG_ONLY)){ // contig and positional range
		const VcfContig* contig = hdr.GetContig(s);
		if(contig == nullptr){
			std::cerr << utility::timestamp("ERROR","INTERVAL") << "Contig does not exist in string " << s << std::endl;
			return false;
		}
		ivecs[contig->idx].push_back(algorithm::Interval<uint32_t,uint32_t>(0,contig->n_bases,0));

		return true;
	} else
		return false;
}


twk_intervals_two::twk_intervals_two() : n_c(0), itree(nullptr){}

twk_intervals_two::twk_intervals_two(const uint32_t n_contigs) :
	n_c(n_contigs),
	itree(new algorithm::IntervalTree< uint32_t, interval_pair_payload >*[n_c])
{
}

twk_intervals_two::~twk_intervals_two(){
	if(itree != nullptr){
		for(uint32_t i = 0; i < n_c; ++i) delete itree[i];
	}
	delete[] itree;
}

bool twk_intervals_two::Build(std::vector<std::string>& strings,
		   const uint32_t n_contigs,
		   const IndexOutput& index,
		   const VcfHeader& hdr)
{
	if(strings.size() == 0) return true;
	for(uint32_t i = 0; i < n_c; ++i) delete itree[i];
	delete[] itree; itree = nullptr;
	ivecs.clear();

	if(n_contigs == 0) return false;
	n_c = n_contigs;
	ivecs.resize(n_c);

	// Reserve memory for 250 intervals for each rid.
	// This is to make sure that pointers of linked reads
	// doesn't become invalid when resizing these vectors.
	for(int i = 0; i < n_c; ++i) ivecs[i].reserve(250);


	// Parse each interval string.
	for(int i = 0; i < strings.size(); ++i){
		if(this->ParseIntervalString(strings[i], hdr) == false)
			return false;
	}

	// Dedupe internals
	this->Dedupe();

	for(uint32_t i = 0; i < ivecs_internal.size(); ++i){
		for(int j = 0; j < ivecs_internal[i].size(); ++j){
			std::vector<IndexEntryOutput*> idx = index.FindOverlap(i,ivecs_internal[i][j].start,ivecs_internal[i][j].stop);
			overlap_blocks.insert(overlap_blocks.end(), idx.begin(), idx.end());
		}
	}

	itree = new algorithm::IntervalTree< uint32_t, interval_pair_payload >*[n_c];
	for(uint32_t i = 0; i < n_c; ++i){
		//std::cerr << "ivecs-before=" << ivecs[i].size() << std::endl;
		itree[i] = new algorithm::IntervalTree< uint32_t, interval_pair_payload >(ivecs[i]);
		//std::cerr << "ivecs-after=" << ivecs[i].size() << std::endl;
	}
	//exit(1);

	if(overlap_blocks.size() == 0){
		std::cerr << utility::timestamp("ERROR","INTERVAL") << "Found no blocks overlapping the provided range(s)..." << std::endl;
		return(false);
	}

	return(true);
}

bool twk_intervals_two::ParseIntervalStrings(std::vector<std::string>& ivals, VcfHeader& hdr){
	if(ivals.size() == 0)
		return true;

	for(int i = 0; i < ivals.size(); ++i){
		if(this->ParseIntervalString(ivals[i], hdr) == false)
			return false;
	}
	return true;
}

bool twk_intervals_two::ParseIntervalString(const std::string& s, const VcfHeader& hdr){
	if(s.size() == 0)
		return false;

	std::vector<std::string> dual_range = utility::split(s, ',');
	if(dual_range.size() > 2){
		std::cerr << utility::timestamp("ERROR","INTERVAL") << "Illegal format: " << s << std::endl;
		return false;
	}

	if(dual_range.size() == 0){
		std::cerr << utility::timestamp("ERROR","INTERAL") << "Illegal format: " << s << std::endl;
		return false;
	}

	// non-anchored string
	if(dual_range.size() == 1){
		if(std::regex_match(s, TWK_REGEX_CONTIG_RANGE)){ // contig only
			std::vector<std::string> rid = utility::split(s, ':');
			std::vector<std::string> pos;
			if(rid.size() == 2){
				pos = utility::split(rid[1], '-');
			} else {
				std::cerr << utility::timestamp("ERROR","INTERVAL") << "Illegal format: " << s << std::endl;
				return false;
			}
			uint32_t ival_from = (uint32_t)std::atof(pos[0].c_str());
			uint32_t ival_to   = (uint32_t)std::atof(pos[1].c_str());
			const VcfContig* contig = hdr.GetContig(rid[0]);
			if(contig == nullptr){
				std::cerr << utility::timestamp("ERROR","INTERVAL") << "Contig does not exist in string " << s << std::endl;
				return false;
			}
			uint32_t ival_rid = contig->idx;
			ivecs[ival_rid].push_back(interval(ival_from,ival_to,interval_pair_payload(-1,0,0)));

			return true;

		} else if(std::regex_match(s, TWK_REGEX_CONTIG_POSITION)){ // contig and position
			std::vector<std::string> rid = utility::split(s, ':');

			uint32_t ival_from = (uint32_t)std::atof(rid[1].c_str());
			uint32_t ival_to   = (uint32_t)std::atof(rid[1].c_str());
			const VcfContig* contig = hdr.GetContig(rid[0]);
			if(contig == nullptr){
				std::cerr << utility::timestamp("ERROR","INTERVAL") << "Contig does not exist in string " << s << std::endl;
				return false;
			}
			uint32_t ival_rid = contig->idx;
			ivecs[ival_rid].push_back(interval(ival_from,ival_to,interval_pair_payload(-1,0,0)));

			return true;
		} else if(std::regex_match(s, TWK_REGEX_CONTIG_ONLY)){ // contig and positional range
			const VcfContig* contig = hdr.GetContig(s);
			if(contig == nullptr){
				std::cerr << utility::timestamp("ERROR","INTERVAL") << "Contig does not exist in string " << s << std::endl;
				return false;
			}
			ivecs[contig->idx].push_back(interval(0,contig->n_bases,interval_pair_payload(-1,0,0)));

			return true;
		} else
			return false;
	} else {
		// is anchored.
		// not possible to dedupe anchored tuples
		uint32_t mateA_offset = 0, mateB_offset = 0;
		int32_t  mateA_rid = -1, mateB_rid = -1;
		interval *mateA = nullptr, *mateB = nullptr;

		for(int i = 0; i < 2; ++i){
			if(std::regex_match(dual_range[i], TWK_REGEX_CONTIG_RANGE)){ // contig only
				std::vector<std::string> rid = utility::split(dual_range[i], ':');
				std::vector<std::string> pos;

				if(rid.size() == 2){
					pos = utility::split(rid[1], '-');
				} else {
					std::cerr << utility::timestamp("ERROR","INTERVAL") << "Illegal format: " << s << std::endl;
					return false;
				}
				uint32_t ival_from = (uint32_t)std::atof(pos[0].c_str());
				uint32_t ival_to   = (uint32_t)std::atof(pos[1].c_str());
				const VcfContig* contig = hdr.GetContig(rid[0]);
				if(contig == nullptr){
					std::cerr << utility::timestamp("ERROR","INTERVAL") << "Contig does not exist in string " << s << std::endl;
					for(int k = 0; k < hdr.contigs_.size(); ++k){
						std::cerr << hdr.contigs_[k].name << std::endl;
					}
					return false;
				}
				uint32_t ival_rid = contig->idx;
				ivecs[ival_rid].push_back(interval(ival_from,ival_to,interval_pair_payload(-1,0,0)));
				if(i == 0){ mateA = &ivecs[ival_rid][ivecs[ival_rid].size()-1]; mateA_offset = ivecs[ival_rid].size() - 1; mateA_rid = ival_rid; }
				else { mateB = &ivecs[ival_rid][ivecs[ival_rid].size()-1]; mateB_offset = ivecs[ival_rid].size() - 1; mateB_rid = ival_rid; }

			} else if(std::regex_match(dual_range[i], TWK_REGEX_CONTIG_POSITION)){ // contig and position
				std::vector<std::string> rid = utility::split(dual_range[i], ':');

				uint32_t ival_from = (uint32_t)std::atof(rid[1].c_str());
				uint32_t ival_to   = (uint32_t)std::atof(rid[1].c_str());
				const VcfContig* contig = hdr.GetContig(rid[0]);
				if(contig == nullptr){
					std::cerr << utility::timestamp("ERROR","INTERVAL") << "Contig does not exist in string " << s << std::endl;
					return false;
				}
				uint32_t ival_rid = contig->idx;
				ivecs[ival_rid].push_back(interval(ival_from,ival_to,interval_pair_payload(-1,0,0)));

				if(i == 0){ mateA = &ivecs[ival_rid][ivecs[ival_rid].size()-1]; mateA_offset = ivecs[ival_rid].size() - 1; mateA_rid = ival_rid; }
				else { mateB = &ivecs[ival_rid][ivecs[ival_rid].size()-1]; mateB_offset = ivecs[ival_rid].size() - 1; mateB_rid = ival_rid; }

			} else if(std::regex_match(dual_range[i], TWK_REGEX_CONTIG_ONLY)){ // contig and positional range
				const VcfContig* contig = hdr.GetContig(dual_range[i]);
				if(contig == nullptr){
					std::cerr << utility::timestamp("ERROR","INTERVAL") << "Contig does not exist in string " << s << std::endl;
					return false;
				}
				ivecs[contig->idx].push_back(interval(0,contig->n_bases,interval_pair_payload(-1,0,0)));
				if(i == 0){ mateA = &ivecs[contig->idx].back(); mateA_offset = ivecs[contig->idx].size() - 1; mateA_rid = contig->idx; }
				else { mateB = &ivecs[contig->idx].back(); mateB_offset = ivecs[contig->idx].size() - 1; mateB_rid = contig->idx; }
			} else {
				return false;
			}
		} // end loop

		mateA->value = interval_pair_payload(mateB_rid, mateB_offset, 0);
		mateB->value = interval_pair_payload(mateA_rid, mateA_offset, 1);
		//std::cerr << "have=" << mateA->start << "-" << mateA->stop << " with " << mateA->value.rid << "," << mateA->value.offset << " and " <<
		//		mateB->start << "-" << mateB->stop << " with " << mateB->value.rid << "," << mateB->value.offset << std::endl;
		return true;
	}
}

void twk_intervals_two::Dedupe(){
	ivecs_internal = ivecs;

	for(int i = 0; i < n_c; ++i){ // every rid
		if(ivecs_internal[i].size() == 0) continue;
		std::vector<interval> ivals;
		std::sort(ivecs_internal[i].begin(), ivecs_internal[i].end());
		ivals.push_back(ivecs_internal[i][0]);
		for(int j = 1; j < ivecs_internal[i].size(); ++j){
			if(  ivecs_internal[i][j].start <= ivals.back().stop
			  && ivecs_internal[i][j].stop  >= ivals.back().start)
			{
				ivals.back().stop = ivecs_internal[i][j].stop;
			}
			else ivals.push_back(ivecs_internal[i][j]);
		}
		ivecs_internal[i] = ivals;
	}


	/*for(int i = 0; i < n_c; ++i){
		for(int j = 0; j < ivecs_internal[i].size(); ++j){
			std::cerr << "ivecs=" << ivecs_internal[i][j].start << "-" << ivecs_internal[i][j].stop << std::endl;
		}
	}*/
}

bool twk_intervals_two::FilterInterval(const twk1_two_t& rec) const {
	std::vector< interval > intervals = this->itree[rec.ridA]->findOverlapping(rec.Apos, rec.Apos);
	uint32_t n_linked = 0, matches_F = 0, matches = 0;

	for(int i = 0; i < intervals.size(); ++i){
		if(intervals[i].value.mate != 0)
			continue;

		++matches_F;
		n_linked += intervals[i].value.rid >= 0;
		if(intervals[i].value.rid >= 0){
			const interval& mate = ivecs[intervals[i].value.rid][intervals[i].value.offset];
			if(rec.Bpos <= mate.stop && rec.Bpos >= mate.start && rec.ridB == intervals[i].value.rid){
				++matches;
				goto finished;
			}
		}
	}

	finished:
	if(n_linked) return(matches == 0);
	return(matches_F == 0);
}

}
