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
	double price = 0;
	/*for (size_t i = 1; i < _paymentSchedule.size(); i++)
	{
		//double dayFractionBetweenFixings = // (Ti)' - (Ti-1)' ... (The (Tj)' are fixing dates)
			//(_fixingSchedule[i].minus(_fixingSchedule[i - 1], _DayCountConvention)).ConvertToDayFraction();
		double dayFraction =
			(_periodsSchedule[i].minus(_periodsSchedule[i - 1], _DayCountConvention)).ConvertToDayFraction();
		double dayFractionPayment = // Ti - t ... (Ti is a payment date and t = 0)
			(_paymentSchedule[i].minus(_startDate, _DayCountConvention)).ConvertToDayFraction();
		double dayFractionPaymentPrevious = // Ti-1 - t ... (Ti-1 is a payment date and t = 0)
			(_paymentSchedule[i - 1].minus(_startDate, _DayCountConvention)).ConvertToDayFraction();
		//double dayFractionPrevious = // (Ti-1)' - t ... (t = 0)
			//(_fixingSchedule[i - 1].minus(_startDate, _DayCountConvention)).ConvertToDayFraction();
		//double dayFractionCurrent = // (Ti)' - t ... (t = 0)
			//(_fixingSchedule[i].minus(_startDate, _DayCountConvention)).ConvertToDayFraction();

		price += ZC(Interpolation(dayFractionPayment, _discount), dayFractionPayment)
			* dayFraction
			* (ZC(Interpolation(dayFractionPayment, _discount), dayFractionPayment)
				/ ZC(Interpolation(dayFractionPaymentPrevious, _discount), dayFractionPaymentPrevious) - 1)
			/ dayFraction;
	}*/
	size_t n = _paymentSchedule.size();
	double d = (_paymentSchedule[n - 1].minus(_startDate, _DayCountConvention)).ConvertToDayFraction();
	price = ZC(Interpolation(0, _discount), 0) - ZC(Interpolation(d, _discount), d);
	return price;
}

