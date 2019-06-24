#ifndef CFIXEDLEG_HH
#define CFIXEDLEG_HH
#include "Functions.hh"
#include "cDate.hh"
#include "cPeriod.hh"
#include "cLeg.hh"

class cFixedLeg : public cLeg
{
private:
	const double _fixedRate;

public:
	cFixedLeg(const cLeg& leg, const double& fixedRate);
	cFixedLeg(const cFixedLeg& fixedLeg);
};









#endif