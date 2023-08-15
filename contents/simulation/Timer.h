#ifndef TIMER
#define TIMER

#include <chrono>


class Timer
{
	/*
	Class to time different functions and code segments.
	*/
private:
	// Type aliases to make accessing nested type easier
	using clock_type = std::chrono::steady_clock;
	using second_type = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<clock_type> m_beg{ clock_type::now() };

public:
	void reset();

	double elapsed() const;
};

#endif
