#ifndef TWK_LD_BALANCING_H_
#define TWK_LD_BALANCING_H_

#include "utility.h"

namespace tomahawk {

/**<
 * Load balancer for calculating linkage-disequilibrium. Partitions the total
 * problem into psuedo-balanced subproblems. The size and number of sub-problems
 * can be parameterized.
 */
struct twk_ld_balancer {
	twk_ld_balancer() : diag(false), n(0), p(0), c(0), fromL(0), toL(0), fromR(0), toR(0), n_m(0){}

	/**<
	 * Find the desired target subproblem range as a tuple (fromL,toL,fromR,toR).
	 * @param n_blocks      Total number of blocks.
	 * @param desired_parts Desired number of subproblems to solve.
	 * @param chosen_part   Target subproblem we are interested in getting the ranges for.
	 * @return              Return TRUE upon success or FALSE otherwise.
	 */
	bool Build(uint32_t n_blocks,
	           uint32_t desired_parts,
	           uint32_t chosen_part)
	{
		if(chosen_part >= desired_parts){
			std::cerr << utility::timestamp("ERROR","BALANCER") << "Illegal chosen block: " << chosen_part << " >= " << desired_parts << std::endl;
			return false;
		}

		n = n_blocks; p = desired_parts; c = chosen_part;
		if(p > n){
			std::cerr << utility::timestamp("ERROR","BALANCER") << "Illegal desired number of blocks! You are asking for more subproblems than there are blocks available (" << p << ">" << n << ")..." << std::endl;
			return false;
		}

		if(p == 1){
			p = 1; c = 0;
			fromL = 0; toL = n_blocks; fromR = 0; toR = n_blocks;
			n_m = n_blocks; diag = true;
			return true;
		}

		uint32_t factor = 0;
		for(uint32_t i = 1; i < desired_parts; ++i){
			if( (((i*i) - i) / 2) + i == desired_parts ){
				//std::cerr << "factor is " << i << std::endl;
				factor = i;
				break;
			}
		}

		if(factor == 0){
			std::cerr << utility::timestamp("ERROR","BALANCER") << "Could not partition into " << desired_parts << " number of subproblems. This number is not a function of x!2 + x..." << std::endl;
			return false;
		}

		// cycle
		uint32_t chunk_size = n / factor;
		uint32_t fL = 0, tL = 0, fR = 0, tR = 0;
		for(uint32_t i = 0, k = 0; i < factor; ++i){ // rows
			for(uint32_t j = i; j < factor; ++j, ++k){ // cols
				tR = (j + 1 == factor ? n_blocks : chunk_size*(j+1));
				fR = tR - chunk_size;
				tL = (i + 1 == factor ? n_blocks : chunk_size*(i+1));
				fL = tL - chunk_size;

				//std::cerr << fL << "-" << tL << "->" << fR << "-" << tR << " total=" << n_blocks << " chunk=" << chunk_size << "desired=" << desired_parts << std::endl;
				if(k == chosen_part){
					//std::cerr << "chosen part:" << std::endl;
					fromL = fL; toL = tL; fromR = fR; toR = tR;
					n_m = (toL - fromL) + (toR - fromR); diag = false;
					if(i == j){ n_m = toL - fromL; diag = true; }
					return true;
				}
			}
		}
		return true;
	}

	/**<
	 * Construction in the special case of (1,...,n) sites compared against
	 * all other sites. This is useful in the case you want to compute a single
	 * site vs all-others very quickly.
	 * @param n_blocks      Total number of blocks.
	 * @param desired_parts Desired number of subproblems to solve.
	 * @param chosen_part   Target subproblem we are interested in getting the ranges for.
	 * @return              Return TRUE upon success or FALSE otherwise.
	 */
	bool BuildSingleSite(uint32_t n_blocks,
	                     uint32_t desired_parts,
	                     uint32_t chosen_part)
	{
		if(desired_parts != 1) return false;
		if(chosen_part > desired_parts) return false;
		p = 1; c = 0;
		fromL = 0; toL = 1; fromR = 0; toR = n_blocks;
		n_m = n_blocks; diag = false;
		return true;
	}

public:
	bool diag; // is selectd chunk diagonal
	uint32_t n, p, c; // number of available blocks, desired parts, chosen part
	uint32_t fromL, toL, fromR, toR;
	uint32_t n_m; // actual blocks used
};

/**<
 * Work balancer for twk_ld_engine threads. Uses a non-blocking spinlock to produce a
 * tuple (from,to) of integers representing the start ref block and dst block.
 * This approach allows perfect load-balancing at a small overall CPU cost. This is
 * directly equivalent to dynamic load-balancing (out-of-order execution) when n > 1.
 */
struct twk_ld_dynamic_balancer {
	typedef bool (twk_ld_dynamic_balancer::*get_func)(uint32_t& from, uint32_t& to, uint8_t& type);

	twk_ld_dynamic_balancer() :
		diag(false), window(false),
		n_perf(0), i(0), j(0),
		fL(0), tL(0), fR(0), tR(0), l_window(0),
		ldd(nullptr),
		_getfunc(&twk_ld_dynamic_balancer::GetBlockPair)
	{}
	~twk_ld_dynamic_balancer(){}

	void operator=(const twk_ld_balancer& balancer){
		diag = balancer.diag;
		fL   = balancer.fromL;
		tL   = balancer.toL;
		fR   = balancer.fromR;
		tR   = balancer.toR;
		i    = balancer.fromL;
		j    = balancer.fromR;
	}

	void SetWindow(const bool window, const uint32_t l_window){
		this->window   = window;
		this->l_window = l_window;
		SetWindow(window);
	}

	/**<
	 * Parameterisation of window mode: should we check if blocks overlap given
	 * some maximum distance between pairs? See function `GetBlockWindow` for
	 * additional details.
	 * @param yes Set or unset window mode.
	 */
	inline void SetWindow(const bool yes = true){
		window = yes;
		_getfunc = (window ? &twk_ld_dynamic_balancer::GetBlockWindow : &twk_ld_dynamic_balancer::GetBlockPair);
	}

	/**<
	 * Indirection using functional pointer to actual function used. This
	 * allows us to use a singular function without writing multiple
	 * versions of downstream functions.
	 * @param from Row position
	 * @param to   Column position
	 * @param type Diagonal (1) or square (0)
	 * @return     Returns TRUE if it is possible to retrieve a new (x,y)-pair or FALSE otherwise.
	 */
	inline bool Get(uint32_t& from, uint32_t& to, uint8_t& type){ return((this->*_getfunc)(from, to, type)); }

	/**<
	 * Retrieves (x,y)-coordinates from the selected load-balancing subproblem.
	 * This variation also checks if the two blocks (x,y) can have any overlapping
	 * regions given some parameterized maximum distance.
	 * Uses a spin-lock to make this function thread-safe.
	 * @param from Row position
	 * @param to   Column position
	 * @param type Diagonal (1) or square (0)
	 * @return     Returns TRUE if it is possible to retrieve a new (x,y)-pair or FALSE otherwise.
	 */
	bool GetBlockWindow(uint32_t& from, uint32_t& to, uint8_t& type){
		spinlock.lock();

		if(j == tR){
			++i; j = (diag ? i : fR); from = i; to = j; type = 1; ++j;
			if(i == tL){ spinlock.unlock(); return false; }
			++n_perf;
			spinlock.unlock();
			return true;
		}
		if(i == tL){ spinlock.unlock(); return false; }

		// First in tgt block - last in ref block
		if(i != j){
			// check if this (x,y) pair have any overlapping intervals.
			if(ldd[j].blk->rcds[0].pos - ldd[i].blk->rcds[ldd[i].n_rec-1].pos > l_window){
				++i; j = (diag ? i : fR); from = i; to = j; type = 1; ++j;
				spinlock.unlock();
				return true;
			}
		}

		type = (i == j); from = i; to = j;
		++j; ++n_perf;

		spinlock.unlock();

		return true;
	}

	/**<
	 * Retrieves (x,y)-coordinates from the selected load-balancing subproblem.
	 * Uses a spin-lock to make this function thread-safe.
	 * @param from Row position
	 * @param to   Column position
	 * @param type Diagonal (1) or square (0)
	 * @return     Returns TRUE if it is possible to retrieve a new (x,y)-pair or FALSE otherwise.
	 */
	bool GetBlockPair(uint32_t& from, uint32_t& to, uint8_t& type){
		spinlock.lock();

		if(j == tR){ // if current position is at the last column
			++i; j = (diag ? i : fR); from = i; to = j; type = 1; ++j;
			// if current position is at the last row
			if(i == tL){ spinlock.unlock(); return false; }
			++n_perf;
			spinlock.unlock();
			return true;
		}
		// if current position is at the last row
		if(i == tL){ spinlock.unlock(); return false; }
		type = (i == j); from = i; to = j;
		++j; ++n_perf;

		spinlock.unlock();

		return true;
	}

public:
	bool diag, window;
	uint32_t n_perf, i,j;
	uint32_t fL, tL, fR, tR, l_window;
	twk1_ldd_blk* ldd;
	get_func _getfunc;
	SpinLock spinlock;
};

}



#endif /* LIB_LD_LD_BALANCING_H_ */
