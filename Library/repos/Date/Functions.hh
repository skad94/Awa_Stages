#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH
#include <vector>
#include <iostream>
#include "cDate.hh"
#include "cPeriod.hh"

using namespace std;

vector<cDate>
Schedule(cDate& start, const cPeriod& maturity, const cPeriod& freq);

void
ShowSchedule(vector<cDate>& schedule);

#endif
