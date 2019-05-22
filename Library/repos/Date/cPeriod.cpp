#include "cPeriod.hh"

cPeriod::cPeriod() : cDate(0, 0, 0), _convention(conv_30_360) {}

cPeriod::cPeriod(const cDate& date, const Convention& convention)
	: cDate(date), _convention(convention) {}

cPeriod::cPeriod(const cPeriod& period)
	: cDate(period), _convention(period._convention) {}


