#include "cDate.hh"

class cPeriod : public cDate
{
private:
	int _monthConvention;
	int _yearConvention;
public:
	cPeriod();
	cPeriod(const int& day, const int& month, const int& year, const int& monthConvention, const int& yearConvention);
	cPeriod(const cPeriod& period);

};
