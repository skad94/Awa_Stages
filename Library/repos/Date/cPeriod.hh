#ifndef CPERIOD_HH
#define CPERIOD_HH

#include "cDate.hh"

enum eConvention { conv_30_360, conv_30_365, conv_ACT_ACT };

class cPeriod : public cDate
{
private:
	eConvention _convention;

public:
	cPeriod();
	cPeriod(const int& day, const int& month, const int& year, const eConvention& convention);
	cPeriod(const cPeriod& period);
	bool IsValid() const override;

};

#endif