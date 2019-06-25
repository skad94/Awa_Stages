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
	for (double i = _freq.ConvertToDayFraction();
		i <= _maturity.ConvertToDayFraction(); i = i + _freq.ConvertToDayFraction())
	{
		price -= ZC(Interpolation(i, _discount), i) * _freq.ConvertToDayFraction() * _fixedRate;
	}
	return price;
}