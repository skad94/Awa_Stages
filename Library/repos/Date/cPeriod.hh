#ifndef CPERIOD_HH
#define CPERIOD_HH

#include "cDate.hh"

class cPeriod : public cDate
{
private:
	int _monthConvention; //eg 30
	int _yearConvention; //eg 360
	//if ACT convention then 0, TODO
public:
	cPeriod();
	cPeriod(const int& day, const int& month, const int& year,
		const int& monthConvention, const int& yearConvention);
	cPeriod(const cPeriod& period);
	int GetMonthConvention() const;
	int GetYearConvention() const;
};

#endif