#include "cPeriod.hh"

cPeriod::cPeriod() : cDate(-1, -1, -1), _convention(conv_30_360) {}

cPeriod::cPeriod(const int& day, const int& month, const int& year, const Convention& convention)
	: cDate(day, month, year), _convention(convention) {}

cPeriod::cPeriod(const cPeriod& period)
	: cDate(period), _convention(period._convention) {}

bool
cPeriod::IsValid() const
{//Checking whether the period is valid 
	return (0 <= _month && _month <= 11 && 0 <= _day && _day <= 30);
}
