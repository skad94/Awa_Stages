#include "cFixedLeg.hh"

cFixedLeg::cFixedLeg(const cLeg& leg, const double& fixedRate) :
	cLeg(leg),
	_fixedRate(fixedRate) {}

cFixedLeg::cFixedLeg(const cFixedLeg& fixedLeg) :
	cLeg(fixedLeg),
	_fixedRate(fixedLeg._fixedRate) {}

double cFixedLeg::PriceLeg() const
{//Price of the fixed leg
	double price = 0;
	for (size_t i = 1; i < _paymentSchedule.size(); i++)
	{
		double dayFraction = 
			(_paymentSchedule[i].minus(_startDate, _DayCountConvention)).ConvertToDayFraction();
		price += ZC(Interpolation(dayFraction, _discount), dayFraction) 
			* _freq.ConvertToDayFraction() * _fixedRate;
	}
	return price;
}