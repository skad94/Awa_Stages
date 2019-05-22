#ifndef CPERIOD_HH
#define CPERIOD_HH

#include "cDate.hh"

enum Convention { conv_30_360, conv_30_365, conv_ACT_ACT };

class cPeriod : public cDate
{
private:
	Convention _convention;

public:
	cPeriod();
	cPeriod(const cDate& date, const Convention& convention);
	cPeriod(const cPeriod& period);
};

#endif