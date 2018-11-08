#ifndef INTERVALS_H_
#define INTERVALS_H_

#include "vcf_utils.h"
#include "intervalTree.h"

namespace tomahawk {

struct interval_pair_payload {
	interval_pair_payload() : rid(-1), offset(0), mate(0){}
	interval_pair_payload(const int32_t rid, const uint32_t mate_offset, const uint8_t type) :
		mate(type), rid(rid), offset(mate_offset)
	{}

	bool operator<(const interval_pair_payload& other) const {
		if(rid < other.rid) return true;
		if(other.rid < rid) return false;
		if(offset < other.offset) return true;
		if(other.offset < offset) return false;
		return true;
	}

	uint8_t mate;
	int32_t rid;
	uint32_t offset;
};

/**<
 * Basic intervals container. Handles tuples of (rid,from_pos,to_pos) and maps
 * them to the local tomahawk index. Have functionality for sorting, deduping,
 * and mapping interval tuples and the resulting index objects.
 */
class twk_intervals {
public:
	// interval (from,to,offset to mate) with rid being implicit in the tree
	typedef algorithm::Interval<uint32_t,uint32_t> interval;

	twk_intervals() : n_c(0), itree(nullptr){}
	twk_intervals(const uint32_t n_contigs) :
		n_c(n_contigs),
		itree(new algorithm::IntervalTree<uint32_t,uint32_t>*[n_c])
	{
		// Reserve memory for 250 intervals for each rid.
		// This is to make sure that pointers of linked reads
		// doesn't become invalid when resizing these vectors.
		for(int i = 0; i < n_c; ++i) ivecs[i].reserve(250);
	}

	virtual ~twk_intervals(){
		for(uint32_t i = 0; i < n_c; ++i) delete itree[i];
		delete[] itree;
	}

	/**<
	 * Dedupes interval vectors tuples (rid,fromA,fromB) by extension. If two
	 * intervals partially overlap in the right end then we extend the left-most
	 * interval to the right-most intervals end position.
	 */
	void Dedupe(){
		for(int i = 0; i < n_c; ++i){ // every rid
			if(ivecs[i].size() == 0) continue;
			std::vector<interval> ivals;
			std::sort(ivecs[i].begin(), ivecs[i].end());
			ivals.push_back(ivecs[i][0]);
			for(int j = 1; j < ivecs[i].size(); ++j){
				if(  ivecs[i][j].start <= ivals.back().stop
				  && ivecs[i][j].stop  >= ivals.back().start)
				{
					ivals.back().stop = ivecs[i][j].stop;
				}
				else ivals.push_back(ivecs[i][j]);
			}
			ivecs[i] = ivals;
		}
	}

public:
	uint32_t n_c;
	std::vector< std::vector< interval > > ivecs; // vector of vectors of intervals
	algorithm::IntervalTree<uint32_t,uint32_t>** itree; // interval tree array
};

class twk_intervals_twk : public twk_intervals {
public:
	twk_intervals_twk() = default;
	twk_intervals_twk(const uint32_t n_contigs) : twk_intervals(n_contigs){}

	/**<
	 * Wrapper function for mapping the internal interval tuples to the index
	 * entries provided by the twk_reader object. If construction fails or
	 * no map hits are found we return FALSE.
	 * @param reader Source twk_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool Build(const uint32_t n_contigs, const Index& index){
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

	/**<
	 * Convenince wrapper for interval strings. Iterate over a vector of unparsed
	 * interval strings and parse them.
	 * @param ivals  Source vector of unparsed interval strings.
	 * @param reader Reference instance of twk_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseIntervalStrings(std::vector<std::string>& ivals, io::VcfHeader& hdr){
		if(ivals.size() == 0)
			return true;

		for(int i = 0; i < ivals.size(); ++i){
			if(this->ParseIntervalString(ivals[i], hdr) == false)
				return false;
		}
		return true;
	}

	/**<
	 * Internal function for parsing a source interval string into a (rid,posA,posB)
	 * tuple. Input string has to match any of the following patterns:
	 *    1) ^rid$
	 *    2) ^rid:pos$
	 *    3) ^rid:pos-pos$
	 *
	 * Numerical values may be expressed in any standard form: for example, 10e6
	 * and 20E6 and 100000. This function will return TRUE if parsing is syntactically
	 * legal and the contig identifier exists in the header. This function will not check
	 * for out-of-bounds errors such as providing an interval exceeding the length of a
	 * contig.
	 * @param s      Source unparsed interval string/
	 * @param reader Reference instance of twk_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseIntervalString(const std::string& s, const io::VcfHeader& hdr){
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
			uint32_t ival_to   = (uint32_t)std::atof(rid[1].c_str());
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

public:
	std::vector<IndexEntry*> overlap_blocks; // overlapping blocks of interest
};

class twk_intervals_two {
public:
	// interval (from,to,offset to mate) with rid being implicit in the tree
	typedef algorithm::Interval< uint32_t, interval_pair_payload > interval;

	twk_intervals_two() : n_c(0), itree(nullptr){}

	twk_intervals_two(const uint32_t n_contigs) :
		n_c(n_contigs),
		itree(new algorithm::IntervalTree< uint32_t, interval_pair_payload >*[n_c])
	{
	}

	~twk_intervals_two(){
		for(uint32_t i = 0; i < n_c; ++i) delete itree[i];
		delete[] itree;
	}

	/**<
	 * Wrapper function for mapping the internal interval tuples to the index
	 * entries provided by the twk_reader object. If construction fails or
	 * no map hits are found we return FALSE.
	 * @param reader Source twk_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool Build(std::vector<std::string>& strings,
	           const uint32_t n_contigs,
	           const IndexOutput& index,
	           const io::VcfHeader& hdr)
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

	/**<
	 * Internal function for parsing a source interval string into a (rid,posA,posB)
	 * tuple. Input string has to match any of the following patterns:
	 *    1) ^rid$
	 *    2) ^rid:pos$
	 *    3) ^rid:pos-pos$
	 *    4) ^rid;rid$
	 *    5) ^rid:pos,rid$
	 *    6) ^rid:pos-pos,rid$
	 *    7) ^rid,rid:pos$
	 *    8) ^rid,rid:pos-pos$
	 *
	 * Numerical values may be expressed in any standard form: for example, 10e6
	 * and 20E6 and 100000. This function will return TRUE if parsing is syntactically
	 * legal and the contig identifier exists in the header. This function will not check
	 * for out-of-bounds errors such as providing an interval exceeding the length of a
	 * contig.
	 * @param s      Source unparsed interval string
	 * @param reader Reference instance of two_reader object.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseIntervalString(const std::string& s, const io::VcfHeader& hdr){
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

	/**<
	 * Dedupes interval vectors tuples (rid,fromA,fromB) by extension. If two
	 * intervals partially overlap in the right end then we extend the left-most
	 * interval to the right-most intervals end position.
	 */
	void Dedupe(){
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

	bool FilterInterval(const twk1_two_t& rec) const {
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

public:
	uint32_t n_c;
	std::vector< std::vector< interval > > ivecs; // vector of vectors of intervals
	std::vector< std::vector< interval > > ivecs_internal; // deduped for finding overlapping indices
	algorithm::IntervalTree< uint32_t, interval_pair_payload >** itree; // interval tree array
	std::vector<IndexEntryOutput*> overlap_blocks; // overlapping blocks of interest
};

}


#endif /* INTERVALS_H_ */
