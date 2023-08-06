#include "Timer.h"


void Timer::reset()
{
	/*
	Reset the timer.
	*/
	m_beg = clock_type::now();
}

double Timer::elapsed() const
{
	/*
	Determine the time elapsed.
	*/
	return std::chrono::duration_cast<second_type>(clock_type::now() - m_beg).count();
}
