#ifndef CFLOATINGLEG_HH
#define CFLOATINGLEG_HH
#include "Functions.hh"
#include "cDate.hh"
#include "cPeriod.hh"

class cFloatingLeg
{
private:
	const eConvention_NonBusinessDay _NonBusinessDayConvention;
	const cPeriod _maturity;
	const cPeriod _freq;
	const cDate _startDate;
	const cPeriod _fixingGap;
	const cPeriod _paymentGap;
	vector<cDate> _fixingSchedule;
	vector<cDate> _paymentSchedule;
	vector<cDate> _periodsSchedule;

public:
	//cFloatingLeg(); //Default constructor for this class makes no sense... ?
	cFloatingLeg(const eConvention_NonBusinessDay& NonBusinessDayConvention, const cPeriod& maturity,
		const cPeriod& freq, const cDate& startDate, const cPeriod& fixingGap, const cPeriod& paymentGap);
	cFloatingLeg(const cFloatingLeg& floatingLeg);

};













#endif
