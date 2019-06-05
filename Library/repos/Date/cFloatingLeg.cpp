#include "cFloatingLeg.hh"

cFloatingLeg::cFloatingLeg(
	const eConvention_NonBusinessDay& NonBusinessDayConvention,
	const cPeriod& maturity,
	const cPeriod& freq,
	const cDate& startDate,
	const cPeriod& fixingGap,
	const cPeriod& paymentGap):
	_NonBusinessDayConvention(NonBusinessDayConvention),
	_maturity(maturity),
	_freq(freq),
	_startDate(startDate),
	_fixingGap(fixingGap),
	_paymentGap(paymentGap)
{
	//Here we have to construct the 3 schedules
	_fixingSchedule = Schedule(_startDate, _maturity, _freq);
	_paymentSchedule = _fixingSchedule;
	_periodsSchedule = _fixingSchedule;
	for (size_t dateIndex = 0; dateIndex < _fixingSchedule.size(); dateIndex++)
	{
		_fixingSchedule[dateIndex] -= _fixingGap;
		_paymentSchedule[dateIndex] += _paymentGap;
	}
	//NOW WE HAVE TO CHECK THE WEEK-ENDS !!!
}

cFloatingLeg::cFloatingLeg(const cFloatingLeg& floatingLeg):
	_NonBusinessDayConvention(floatingLeg._NonBusinessDayConvention),
	_maturity(floatingLeg._maturity),
	_freq(floatingLeg._freq),
	_startDate(floatingLeg._startDate),
	_fixingGap(floatingLeg._fixingGap),
	_paymentGap(floatingLeg._paymentGap),
	_fixingSchedule(floatingLeg._fixingSchedule),
	_paymentSchedule(floatingLeg._paymentSchedule),
	_periodsSchedule(floatingLeg._paymentSchedule) {}