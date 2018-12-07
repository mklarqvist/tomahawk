#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <fstream>

namespace tomahawk {

/**<
 * Simple timer class for tracking time differences between two timepoints.
 * Internally use the `chrono::high_resolution_clock` struct such that we can
 * track very short time frames.
 */
class Timer {
public:
	explicit Timer(){}

	/**<
	 * Start the timer by setting the current timestamp as the reference
	 * time.
	 */
	void Start(void){ this->_start = std::chrono::high_resolution_clock::now(); }

	/**<
	 * Returns the number `chrono::duration<double>` object for time elapsed.
	 * If you are interested in the number of seconds elapsed then chain this
	 * function with the child function `count`: `timer.Elapsed().count()`.
	 * @return
	 */
	inline std::chrono::duration<double> Elapsed() const{
		return(std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - this->_start));
	}

	friend std::ostream& operator<<(std::ostream& out, const Timer& timer){
		return out << timer.Elapsed().count();
	}

	std::string ElapsedString(void){ return this->SecondsToTimestring(this->Elapsed().count()); }

private:
	std::string SecondsToTimestring(const double seconds){
		const int32_t hours     = ((int32_t)seconds / 60 / 60);
		const int32_t minutes   = ((int32_t)seconds / 60) % 60;
		const int32_t sec       = (int32_t)seconds % 60;
		const int32_t remainder = (seconds - (int32_t)seconds)*1000;

		if(hours > 0){
			sprintf(&this->buffer[0], "%02uh%02um%02u,%03us",
					hours,
					minutes,
					sec,
					remainder);

			return(std::string(&this->buffer[0], 13));
		} else if(minutes > 0){
			sprintf(&this->buffer[0], "%02um%02u,%03us",
					minutes,
					sec,
					remainder);

			return(std::string(&this->buffer[0], 10));
		} else {
			sprintf(&this->buffer[0], "%02u,%03us",
					sec,
					remainder);

			return(std::string(&this->buffer[0], 7));
		}
	}

private:
	char buffer[64];
	std::chrono::high_resolution_clock::time_point _start;

};

}
#endif /* TIMER_H_ */
