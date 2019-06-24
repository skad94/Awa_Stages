#include "cFixedLeg.hh"

cFixedLeg::cFixedLeg(const cLeg& leg, const double& fixedRate) :
	cLeg(leg),
	_fixedRate(fixedRate) {}

cFixedLeg::cFixedLeg(const cFixedLeg& fixedLeg) :
	cLeg(fixedLeg),
	_fixedRate(fixedLeg._fixedRate) {}