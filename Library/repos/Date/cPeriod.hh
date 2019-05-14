#include "cDate.hh"

class cPeriod : public cDate
{
private:
	int _monthConvention;
	int _yearConvention;
public:
	cPeriod();
	cPeriod(int day, int month, int year, int monthConvention, int yearConvention);
	cPeriod(const cPeriod& period);
	~cPeriod();

};
