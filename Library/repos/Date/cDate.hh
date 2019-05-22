#ifndef CDATE_HH
#define CDATE_HH
#include <iostream>
#include <vector>
#include <string>
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
		void ShowDate();
		bool IsValid() const;
		void SetAsInvalid();
		static bool IsLeapYear(const int& year);
		cDate& operator-=(const cPeriod& period);
		cDate operator-(const cPeriod& period) const;
		cDate& operator+=(const cPeriod& period);
		cDate operator+(const cPeriod& period) const;
		bool operator<(const cDate& d1);
};

#endif



