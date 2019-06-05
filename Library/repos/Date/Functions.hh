#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH
#include <vector>
#include <iostream>
#include "cDate.hh"
#include "cPeriod.hh"

using namespace std;

vector<cDate>
Schedule(const cDate& start, const cPeriod& maturity, const cPeriod& freq);

void
ShowSchedule(const vector<cDate>& schedule);

cDate NumberOfDays_To_Date(int ndays);

int Date_To_NumberOfDays(const cDate& date);

#endif
