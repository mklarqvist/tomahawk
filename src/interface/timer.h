#ifndef TIMER_H_
#define TIMER_H_

#include <string>
#include <iostream>
#include <chrono>

namespace tomahawk{
namespace interface{

class Timer {
public:
	explicit Timer(){}

	void Start(void){ this->_start = std::chrono::high_resolution_clock::now(); }

	std::chrono::duration<double> Elapsed() const{
		return std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - this->_start);
	}

	template <typename T, typename Traits>
	friend std::basic_ostream<T, Traits>& operator<<(std::basic_ostream<T, Traits>& out, const Timer& timer){
		return out << timer.Elapsed().count();
	}

	std::string ElapsedString(void){ return this->SecondsToTimestring(this->Elapsed().count()); }

private:
	std::string SecondsToTimestring(const double seconds){
		const S32 hours = ((S32)seconds / 60 / 60);
		const S32 minutes = ((S32)seconds / 60) % 60;
		const S32 sec = (S32)seconds % 60;
		const S32 remainder = (seconds - (S32)seconds)*1000;

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
	char buffer[13];
	std::chrono::high_resolution_clock::time_point _start;

};

}
}
#endif /* TIMER_H_ */
