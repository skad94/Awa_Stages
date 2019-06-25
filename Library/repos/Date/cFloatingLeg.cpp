#include "cFloatingLeg.hh"

cFloatingLeg::cFloatingLeg(const cLeg& leg, const cPeriod& fixingGap):
	cLeg(leg),
	_fixingGap(fixingGap)
{
	_fixingSchedule = _periodsSchedule;
	for (size_t dateIndex = 0; dateIndex < _periodsSchedule.size(); dateIndex++) 
	{
		_fixingSchedule[dateIndex] -= _fixingGap;
	}
	SetAsValidSchedule(_fixingSchedule, _NonBusinessDayConvention);
}

cFloatingLeg::cFloatingLeg(const cFloatingLeg& floatingLeg):
	cLeg(floatingLeg),
	_fixingGap(floatingLeg._fixingGap),
	_fixingSchedule(floatingLeg._fixingSchedule) {}

double cFloatingLeg::PriceLeg() const
{//Price of the floating leg
	return ZC(Interpolation(0, _discount), 0) 
		- ZC(Interpolation(_maturity.ConvertToDayFraction(), _discount), _maturity.ConvertToDayFraction());
}
