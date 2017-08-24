#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>

namespace Tomahawk{
namespace Interface{

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

	std::string ElapsedPretty(void) const{ return this->SecondsToTimestring(this->Elapsed().count()); }

private:
	std::string SecondsToTimestring(const double seconds) const{
		const S32 hours = ((S32)seconds / 60 / 60);
		const S32 minutes = ((S32)seconds / 60) % 60;
		const S32 sec = (S32)seconds % 60;
		const S32 remainder = (seconds - (S32)seconds)*1000;

		char pad = 0, pad2 = 0;
		if(remainder < 10){ pad = '0'; pad2 = '0';}
		else if(remainder < 100){pad = '0';}

		std::stringstream st;

		if(hours > 0) st << hours << 'h';
		if(minutes > 0) st << minutes << 'm';
		st << sec << '.' << pad << pad2 << remainder << 's';

		return st.str();
	}

private:
	std::chrono::high_resolution_clock::time_point _start;

};

}
}
#endif /* TIMER_H_ */
