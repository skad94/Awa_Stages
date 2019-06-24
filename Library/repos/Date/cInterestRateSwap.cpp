#include "cInterestRateSwap.hh"

cInterestRateSwap::cInterestRateSwap(const cLeg& firstLeg, const cLeg& secondLeg) :
	_firstLeg(firstLeg),
	_secondLeg(secondLeg) {}

cInterestRateSwap::cInterestRateSwap(const cInterestRateSwap& interestRateSwap) :
	_firstLeg(interestRateSwap._firstLeg),
	_secondLeg(interestRateSwap._secondLeg) {}

double cInterestRateSwap::Price_IRS() const
{
	return 0; // Comment Pricer ??
}
