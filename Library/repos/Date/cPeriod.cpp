#include "cPeriod.hh"

cPeriod::cPeriod() :cDate(0, 0, 0), _monthConvention(30), _yearConvention(360) {}

cPeriod::cPeriod(const int& day, 
	             const int& month,
	             const int& year,
	             const int& monthConvention,
	             const int& yearConvention)
	: cDate(day, month, year), _monthConvention(monthConvention), 
	  _yearConvention(yearConvention) {}

cPeriod::cPeriod(const cPeriod& period)
	: cDate(period), _monthConvention(period._monthConvention), _yearConvention(period._yearConvention) {}
