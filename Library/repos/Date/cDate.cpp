#include "cDate.hh"
#include "cPeriod.hh"
#include <iostream>
using namespace std;

cDate::cDate() 
	:_day(1), _month(1), _year(1900) {}

cDate::cDate(const int& day, const int& month, const int& year)
	: _day(day), _month(month), _year(year) {}

cDate::cDate(const cDate& d) 
	: _day(d._day), _month(d._month), _year(d._year) {}

cDate::~cDate() {}

void 
cDate::ShowDate()
{
	cout << _day << "/" << _month << "/" << _year << endl;
}




