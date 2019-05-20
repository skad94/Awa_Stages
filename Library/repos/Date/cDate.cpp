#include "cDate.hh"
#include "cPeriod.hh"

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

bool 
cDate::LeapYear(const int& year)
{//1 if this is a leap year, else 0
	return (year % 400 == 0 || (year % 4 == 0 && year % 100 != 0));
}

cDate& 
cDate::operator-=(const cPeriod& period)
{
	_year -= period._year;
	int monthDiff = _month - period._month;
	if (monthDiff <= 0)
	{
		_year -= 1;
		_month = monthDiff + 12;
	}
	else 
		_month = monthDiff;
	int dayDiff = _day - period._day;
	if (dayDiff <= 0)
	{
		_month -= 1;
		if (_month == 4 || _month == 6 || _month == 9 || _month == 11)
			_day = dayDiff + 30;
		else if (_month == 2)
			if (LeapYear(_year))
				_day = dayDiff + 29;
			else
				_day = dayDiff + 28;
		else
			_day = dayDiff + 30;
	}
	else
		_day = dayDiff;
	return *this;
}

cDate 
cDate::operator-(const cPeriod& period)
{
	cDate res(*this);
	res -= period;
	return res;
}

cDate&
cDate::operator+=(const cPeriod& period)
{
	_year += period._year;
	int monthAdd = _month + period._month;
	if (monthAdd > 12)
	{
		_year += 1;
		_month = monthAdd - 12;
	}
	else
		_month = monthAdd;
	int dayAdd = _day + period._day;
	if ((_month == 4 || _month == 6 || _month == 9 || _month == 11) && dayAdd > 30)
	{
		_month += 1;
		_day = dayAdd - 30;
	}
	else if (_month == 2 && (LeapYear(_year) && dayAdd > 29))
	{
		_month += 1;
		_day = dayAdd - 29;
	}
	else if (_month == 2 && (!LeapYear(_year) && dayAdd > 28))
	{
		_month += 1;
		_day = dayAdd - 28;
	}
	else if (dayAdd > 31)
	{
		_month += 1;
		_day = dayAdd - 31;
	}
	else
		_day = dayAdd;
	return *this;
}

cDate
cDate::operator+(const cPeriod& period)
{
	cDate res(*this);
	res += period;
	return res;
}
