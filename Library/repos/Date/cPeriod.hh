#ifndef CPERIOD_HH
#define CPERIOD_HH

#include "cDate.hh"

class cPeriod : public cDate
{
private:
	eConvention _convention; //Day count convention

public:
	cPeriod();
	cPeriod(const int& day, const int& month, const int& year, const eConvention& convention);
	cPeriod(const cPeriod& period);
	bool IsValid() const override;
	double ConvertToDayFraction() const;

};

#endif