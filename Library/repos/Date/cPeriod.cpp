#include "cPeriod.hh"

cPeriod::cPeriod() : cDate(-1, -1, -1), _convention(conv_30_360) {}

cPeriod::cPeriod(const int& day, const int& month, const int& year, const eConvention& convention)
	: cDate(day, month, year), _convention(convention) {}

cPeriod::cPeriod(const cPeriod& period)
	: cDate(period), _convention(period._convention) {}

bool
cPeriod::IsValid() const
{//Checking whether the period is valid 
	return (0 <= _month && _month <= 11 && 0 <= _day && _day <= 30);
}

double 
cPeriod::ConvertToDayFraction() const
{//Convert a period to a fraction of day according to the day count convention
	if (_convention == conv_30_360)
		return (360 * _year + 30 * _month + _day) / 360.0;
	else if (_convention == conv_30_365)
		return (365 * _year + 30 * _month + _day) / 365.0;
	else
		return 0; //FIX THIS !!!
}