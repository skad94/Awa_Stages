#ifndef CINTERESTRATESWAP_HH
#define CINTERESTRATESWAP_HH
#include "cFloatingLeg.hh"
#include "cFixedLeg.hh"

class cInterestRateSwap 
{
private:
	const cLeg& _firstLeg;
	const cLeg& _secondLeg;

public:
	cInterestRateSwap(const cLeg& firstLeg, const cLeg& secondLeg);
	cInterestRateSwap(const cInterestRateSwap& interestRateSwap);
	double Price_IRS() const;
};








#endif