#ifndef CDATE_HH
#define CDATE_HH
#include <iostream>
using namespace std;

class cPeriod;

class cDate 
{
	private:
		int _day;
		int _month;
		int _year;
	public:
		cDate();
		cDate(const int& day, const int& month, const int& year);
		cDate(const cDate& d);
		~cDate();
		void ShowDate();
		static bool LeapYear(const int& year);
		cDate& operator-=(const cPeriod& period);
		cDate operator-(const cPeriod& period);
		cDate& operator+=(const cPeriod& period);
		cDate operator+(const cPeriod& period);
};

#endif



