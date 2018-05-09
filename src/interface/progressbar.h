#ifndef PROGRESSBAR_H_
#define PROGRESSBAR_H_

#include <iostream>
#include <string>
#include <chrono>
#include <ctime>
#include <sstream>
#include <thread>

#include "../support/helpers.h"
#include "timer.h"

namespace Tomahawk{
namespace Interface{

class ProgressBar {
	typedef ProgressBar self_type;
	typedef Timer timer_type;

public:
	ProgressBar(const bool detailed = false) :
		samples(1),
		totalComparisons(0),
		counter(0),
		outputCount(0),
		Run(true),
		Permission(true),
		Detailed(detailed)
	{}
	~ProgressBar(){}

	void Start(void){
		this->timer.Start();
		this->Run = true;
		this->WorkerThread__ = std::thread(&self_type::__run, this);
		this->WorkerThread__.detach(); // detach tread (do not wait for join)
	}

	void Stop(void){
		this->Run = false;
	}

	void SetDetailed(const bool detailed){ this->Detailed = detailed; }
	void SetSamples(const U64 samples){ this->samples = samples; }
	void SetComparisons(const U64 comparisons){ this->totalComparisons = comparisons; }

	void __run(void){
		U32 i = 0;

		char support_buffer[256];
		std::cerr << Helpers::timestamp("PROGRESS")
				<< std::setw(12) << "Time elapsed"
				<< std::setw(15) << "Variants"
				<< std::setw(20) << "Genotypes"
				<< std::setw(15) << "Output"
				<< std::setw(10) << "Progress"
				<< "\tEst. Time left" << std::endl;
		while(this->Run){
			// Triggered every cycle (119 ms)
			if(this->Detailed)
				this->GetElapsedTime();

			++i;
			if(i % 252 == 0){ // Approximately every 30 sec (30e3 / 119 = 252)
				const double ComparisonsPerSecond = (double)this->counter/this->timer.Elapsed().count();

				const U32 n_p = sprintf(&support_buffer[0], "%0.3f", (double)this->counter.load()/this->totalComparisons*100);
				support_buffer[n_p] = '%';
				std::cerr << Helpers::timestamp("PROGRESS")
						<< std::setw(12) << this->timer.ElapsedString()
						<< std::setw(15) << Helpers::ToPrettyString(this->counter)
						<< std::setw(20) << Helpers::ToPrettyString(this->counter*this->samples)
						<< std::setw(15) << Helpers::ToPrettyString(this->outputCount)
						<< std::setw(10) << std::string(&support_buffer[0], n_p + 1) << '\t'
						<< Helpers::secondsToTimestring((this->totalComparisons - this->counter.load())/ComparisonsPerSecond) << std::endl;
				i = 0;
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(119));
		}

		// Final output
		if(this->Detailed)
			this->GetElapsedTime();

		std::cerr << Helpers::timestamp("PROGRESS") << this->timer.ElapsedString() << "\t" << Helpers::ToPrettyString(this->counter) << "\t" << Helpers::ToPrettyString(this->counter*this->samples) << "\t" << Helpers::ToPrettyString(this->outputCount) << std::endl;
		std::cerr << Helpers::timestamp("PROGRESS") << "Finished" << std::endl;
	}

	void GetElapsedTime(void){
		if(!this->Permission || !this->Run)
			return;

		// Convert to string
		const double ComparisonsPerSecond = (double)this->counter/this->timer.Elapsed().count();
		const std::string ComparisonsPerSecondString = Helpers::NumberThousandsSeparator(std::to_string((U32)ComparisonsPerSecond));

		// Print output
		std::cerr << "\33[2K\r" << this->timer.ElapsedString() << " | Comparisons: " <<
				Helpers::NumberThousandsSeparator(std::to_string(this->counter)) <<
				" (" << ComparisonsPerSecondString << " comparisons/s)" << std::flush;

		return;
	}

	inline const std::atomic<U64>& GetCounter(void) const{ return this->counter; }
	inline const std::atomic<U64>& GetOutputCounter(void) const{ return this->outputCount; }

	inline void operator()(const U64 addCounter, const U32 addOutput){
		this->counter     += addCounter;
		this->outputCount += addOutput;
	}

private:
	U64 samples;
	U64 totalComparisons;
	std::atomic<U64> counter;
	std::atomic<U64> outputCount;
	timer_type timer;
	bool Run;
	bool Permission;
	bool Detailed;
	std::thread WorkerThread__;
};

}
}

#endif /* PROGRESSBAR_H_ */
