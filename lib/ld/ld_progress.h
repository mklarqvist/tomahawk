#ifndef LIB_LD_LD_PROGRESS_H_
#define LIB_LD_LD_PROGRESS_H_

#include <thread>

#include "utility.h"
#include "timer.h"

namespace tomahawk {

/**<
 * Progress ticker for calculating linkage-disequilbirium. Spawns and detaches
 * a thread to tick occasionally in the background. Slaves computing linkage-
 * disequilibrium send their progress to this ticker that collates and summarize
 * that data.
 */
struct twk_ld_progress {
	twk_ld_progress() :
		is_ticking(false), n_s(0), n_cmps(0), n_var(0),
		n_pair(0), n_out(0), b_out(0), thread(nullptr)
	{}
	~twk_ld_progress() = default;

	/**<
	 * Starts the progress ticker. Spawns a detached thread ticking every 30 seconds
	 * in the background until the flag `is_ticking` is set to FALSE or the program
	 * finishes.
	 * @return Returns a pointer to the detached thread.
	 */
	std::thread* Start(){
		delete thread;
		is_ticking = true;
		thread = new std::thread(&twk_ld_progress::StartTicking, this);
		thread->detach();
		return(thread);
	}

	/**<
	 * Internal function displaying the progress message every 30 seconds. This
	 * function is called exclusively by the detached thread.
	 */
	void StartTicking(){
		timer.Start();

		uint64_t variant_overflow = 99E9, genotype_overflow = 999E12;
		uint8_t  variant_width = 15, genotype_width = 20;

		//char support_buffer[256];
		std::cerr << utility::timestamp("PROGRESS")
				<< std::setw(12) << "Time elapsed"
				<< std::setw(variant_width) << "Variants"
				<< std::setw(genotype_width) << "Genotypes"
				<< std::setw(15) << "Output"
				<< std::setw(10) << "Progress"
				<< "\tEst. Time left" << std::endl;

		std::this_thread::sleep_for(std::chrono::seconds(30)); // first sleep
		while(is_ticking){
			if(n_var.load() > variant_overflow)     { variant_width  += 3; variant_overflow  *= 1e3; }
			if(n_var.load()*n_s > genotype_overflow){ genotype_width += 3; genotype_overflow *= 1e3; }

			if(n_cmps){
				std::cerr << utility::timestamp("PROGRESS")
						<< std::setw(12) << timer.ElapsedString()
						<< std::setw(variant_width) << utility::ToPrettyString(n_var.load())
						<< std::setw(genotype_width) << utility::ToPrettyString(n_var.load()*n_s)
						<< std::setw(15) << utility::ToPrettyString(n_out.load())
						<< std::setw(10) << (double)n_var.load()/n_cmps*100 << "%\t"
						<< utility::SecondsToTimestring((n_cmps - n_var.load()) / ((double)n_var.load()/timer.Elapsed().count())) << std::endl;
			} else {
				std::cerr << utility::timestamp("PROGRESS")
						<< std::setw(12) << timer.ElapsedString()
						<< std::setw(variant_width) << utility::ToPrettyString(n_var.load())
						<< std::setw(genotype_width) << utility::ToPrettyString(n_var.load()*n_s)
						<< std::setw(15) << utility::ToPrettyString(n_out.load())
						<< std::setw(10) << 0 << '\t'
						<< 0 << std::endl;
			}
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
		std::cerr << utility::timestamp("PROGRESS") << "Finished in " << this->timer.ElapsedString()
				<< ". Variants: " << utility::ToPrettyString(n_var.load()) << ", genotypes: "
				<< utility::ToPrettyString(n_var.load()*n_s) << ", output: "
				<< utility::ToPrettyString(n_out.load()) << std::endl;
		std::cerr << utility::timestamp("PROGRESS") << utility::ToPrettyString((uint64_t)((double)n_var.load()/timer.Elapsed().count())) << " variants/s and "
				<< utility::ToPrettyString((uint64_t)(((double)n_var.load()*n_s)/timer.Elapsed().count())) << " genotypes/s" << std::endl;
	}

public:
	bool is_ticking;
	uint32_t n_s; // number of samples
	uint64_t n_cmps; // number of comparisons we estimate to perform
	std::atomic<uint64_t> n_var, n_pair, n_out, b_out; // counters used by ld threads
	std::thread* thread; // detached thread
	Timer timer; // timer instance
};

}



#endif /* LIB_LD_LD_PROGRESS_H_ */
