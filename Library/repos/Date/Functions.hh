#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH
#include <vector>
#include <iostream>
//#include "cDate.hh"
//#include "cPeriod.hh"
//#include "cFloatingLeg.hh"

using namespace std;

class cDate;
class cPeriod;
class cFloatingLeg;

enum eConvention_NonBusinessDay { GoForward, GoBackward, GoToTheClosest };

vector<cDate>
Schedule(const cDate& start, const cPeriod& maturity, const cPeriod& freq, const eConvention_NonBusinessDay& NonBusinessDayConvention);

void
ShowSchedule(const vector<cDate>& schedule);

cDate NumberOfDays_To_Date(int ndays);

int Date_To_NumberOfDays(const cDate& date);

#endif
