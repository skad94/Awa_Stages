#ifndef CFLOATINGLEG_HH
#define CFLOATINGLEG_HH
#include "Functions.hh"
#include "cDate.hh"
#include "cPeriod.hh"
#include "cLeg.hh"

class cFloatingLeg : public cLeg
{
private:
	const cPeriod _fixingGap;
	vector<cDate> _fixingSchedule;

public:
	cFloatingLeg(const cLeg& leg, const cPeriod& fixingGap);
	cFloatingLeg(const cFloatingLeg& floatingLeg);
	double PriceLeg() const override;
};













#endif
