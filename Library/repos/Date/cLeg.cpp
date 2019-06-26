#include "cLeg.hh"

cLeg::cLeg(
	const eConvention_NonBusinessDay& NonBusinessDayConvention,
	const eConvention& DayCountConvention,
	const cPeriod& maturity,
	const cPeriod& freq,
	const cDate& startDate,
	const cPeriod& paymentGap,
	const map<double, double> discount):
	_NonBusinessDayConvention(NonBusinessDayConvention),
	_DayCountConvention(DayCountConvention),
	_maturity(maturity),
	_freq(freq),
	_startDate(startDate),
	_paymentGap(paymentGap),
	_discount(discount)
{
	//Here we have to construct the 2 schedules
	_periodsSchedule = Schedule(_startDate, _maturity, _freq, _NonBusinessDayConvention);
	_paymentSchedule = _periodsSchedule;
	for (size_t dateIndex = 0; dateIndex < _periodsSchedule.size(); dateIndex++)
	{
		_paymentSchedule[dateIndex] += _paymentGap;
	}
	SetAsValidSchedule(_periodsSchedule, NonBusinessDayConvention);
	SetAsValidSchedule(_paymentSchedule, NonBusinessDayConvention);
}

cLeg::cLeg(const cLeg& leg):
	_NonBusinessDayConvention(leg._NonBusinessDayConvention),
	_DayCountConvention(leg._DayCountConvention),
	_maturity(leg._maturity),
	_freq(leg._freq),
	_startDate(leg._startDate),
	_paymentGap(leg._paymentGap),
	_periodsSchedule(leg._periodsSchedule),
	_paymentSchedule(leg._paymentSchedule),
	_discount(leg._discount) {}
