#ifndef LIB_SORT_PROGRESS_H_
#define LIB_SORT_PROGRESS_H_

#include <thread>
#include <atomic>

#include "utility.h"
#include "timer.h"

namespace tomahawk {

/**<
 * Progress ticker for calculating linkage-disequilbirium. Spawns and detaches
 * a thread to tick occasionally in the background. Slaves computing linkage-
 * disequilibrium send their progress to this ticker that collates and summarize
 * that data.
 */
struct twk_sort_progress {
	twk_sort_progress() :
		is_ticking(false), n_cmps(0), cmps(0), thread(nullptr)
	{}
	~twk_sort_progress() = default;

	/**<
	 * Starts the progress ticker. Spawns a detached thread ticking every 30 seconds
	 * in the background until the flag `is_ticking` is set to FALSE or the program
	 * finishes.
	 * @return Returns a pointer to the detached thread.
	 */
	std::thread* Start(){
		delete thread;
		is_ticking = true;
		thread = new std::thread(&twk_sort_progress::StartMerge, this);
		thread->detach();
		return(thread);
	}

	/**<
	 * Internal function displaying the progress message every 30 seconds. This
	 * function is called exclusively by the detached thread.
	 */
	void StartMerge(){
		timer.Start();

		uint64_t variant_overflow = 99E9;
		uint8_t  variant_width = 15;

		//char support_buffer[256];
		std::cerr << utility::timestamp("PROGRESS")
				<< std::setw(12) << "Time elapsed"
				<< std::setw(variant_width) << "Variants"
				<< std::setw(10) << "Progress"
				<< "\tEst. Time left" << std::endl;

		std::this_thread::sleep_for(std::chrono::seconds(30)); // first sleep
		while(is_ticking){
			if(cmps.load() > variant_overflow) { variant_width  += 3; variant_overflow  *= 1e3; }

			std::cerr << utility::timestamp("PROGRESS")
					<< std::setw(12) << timer.ElapsedString()
					<< std::setw(variant_width) << utility::ToPrettyString(cmps.load())
					<< std::setw(10) << (double)cmps.load()/n_cmps*100 << "%\t"
					<< utility::SecondsToTimestring((n_cmps - cmps.load()) / ((double)cmps.load()/timer.Elapsed().count())) << std::endl;

			std::this_thread::sleep_for(std::chrono::seconds(30));
		}
	}

	/**<
	 * Print out the final tally of time elapsed, number of variants computed,
	 * and average throughput. This method cannot be made const as the function
	 * ElapsedString in the Timer class internally updates a buffer for performance
	 * reasons. This has no consequence as this function is ever only called once.
	 */
	void PrintFinal(){
		std::cerr << utility::timestamp("PROGRESS") << this->timer.ElapsedString() << "\t"
				<< utility::ToPrettyString(cmps.load())
		        << " (" << utility::ToPrettyString((uint64_t)((double)cmps.load()/timer.Elapsed().count())) << " variants/s)" << std::endl;
		std::cerr << utility::timestamp("PROGRESS") << "Finished!" << std::endl;
	}

public:
	bool is_ticking;
	uint64_t n_cmps; // number of comparisons we estimate to perform
	std::atomic<uint64_t> cmps;
	std::thread* thread; // detached thread
	Timer timer; // timer instance
};

}


#endif /* LIB_SORT_PROGRESS_H_ */
