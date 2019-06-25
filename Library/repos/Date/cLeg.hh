#ifndef CLEG_HH
#define CLEG_HH
#include "Functions.hh"
#include "cDate.hh"
#include "cPeriod.hh"

class cLeg 
{
protected:
	const eConvention_NonBusinessDay _NonBusinessDayConvention;
	const cPeriod _maturity;
	const cPeriod _freq;
	const cDate _startDate;
	const cPeriod _paymentGap;
	vector<cDate> _paymentSchedule;
	vector<cDate> _periodsSchedule;
	const map<double, double> _discount;

public:
	cLeg(const eConvention_NonBusinessDay& NonBusinessDayConvention, const cPeriod& maturity,
		const cPeriod& freq, const cDate& startDate, 
		const cPeriod& paymentGap, const map<double, double> discount);
	cLeg(const cLeg& leg);
	virtual double PriceLeg() const = 0;
};












#endif
