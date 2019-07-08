#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
//#include "cDate.hh"
//#include "cPeriod.hh"
//#include "cFloatingLeg.hh"

using namespace std;

class cDate;
class cPeriod;
class cLeg;
class cFloatingLeg;
class cFixedLeg;
class cInterestRateSwap;

enum eConvention { conv_30_360, conv_30_365, conv_ACT_ACT };
enum eConvention_NonBusinessDay { GoForward, GoBackward, GoToTheClosest };

vector<cDate>
Schedule(const cDate& start, const cPeriod& maturity, const cPeriod& freq, const eConvention_NonBusinessDay& NonBusinessDayConvention);

void
ShowSchedule(const vector<cDate>& schedule);

cDate NumberOfDays_To_Date(int ndays);

int Date_To_NumberOfDays(const cDate& date);

void
SetAsValidSchedule(vector<cDate>& schedule, const eConvention_NonBusinessDay& NonBusinessDayConvention);

double
Interpolation(const double& x, const map<double, double> curve, const string convention = "linear");

double
ZC(const double& rate, const double& maturity, string convention = "composed", const double& t = 0);

void ReplaceComa(string& s);

double ConversionReuters(const string& maturity, const eConvention& convention);

map<double, double> YieldCurve(string name, const eConvention& convention, string user = "alepeltier");

void Affiche(const map<double, double> mymap);

#endif
