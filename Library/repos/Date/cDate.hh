#ifndef CDATE_HH
#define CDATE_HH
#include <iostream>
#include <vector>
#include <string>
#include "Functions.hh"
using namespace std;

class cPeriod;

class cDate 
{
	protected:
		int _day;
		int _month;
		int _year;
	public:
		cDate();
		cDate(const int& day, const int& month, const int& year);
		cDate(const cDate& d);
		int GetDay() const;
		int GetMonth() const;
		int GetYear() const;
		void Show() const;
		virtual bool IsValid() const;
		void SetAsInvalid();
		static bool IsLeapYear(const int& year);
		cDate& operator-=(const cPeriod& period);
		cDate operator-(const cPeriod& period) const;
		cDate& operator+=(const cPeriod& period);
		cDate operator+(const cPeriod& period) const;
		bool operator<(const cDate& d1) const;
		int WhatDayIsIt() const;

};

#endif



