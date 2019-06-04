#include "cDate.hh"
#include "cPeriod.hh"

cDate::cDate() 
	:_day(-1), _month(-1), _year(-1) {}

cDate::cDate(const int& day, const int& month, const int& year)
	: _day(day), _month(month), _year(year) {}

cDate::cDate(const cDate& d) 
	: _day(d._day), _month(d._month), _year(d._year) {}

int 
cDate::GetDay() const
{
	return _day;
}

int
cDate::GetMonth() const
{
	return _month;
}

int
cDate::GetYear() const
{
	return _year;
}

void 
cDate::Show()
{//Display date as day/month/year
	cout << _day << "/" << _month << "/" << _year << endl;
}

bool 
cDate::IsValid() const
{//Checking whether the date is valid 
	if (1 <= _month && _month <= 12)
	{
		if (_month == 4 || _month == 6 || _month == 9 || _month == 11) //Months with 30 days
			return (1 <= _day && _day <= 30);
		else if (_month == 2) //February
			if (IsLeapYear(_year))
				return (1 <= _day && _day <= 29);
			else
				return (1 <= _day && _day <= 28);
		else //Months with 31 days
			return (1 <= _day && _day <= 31);
	}
	return 0;
}

void 
cDate::SetAsInvalid()
{//Set the date as invalid
	_day = -1;
	_month = -1;
	_year = -1;
}

bool 
cDate::IsLeapYear(const int& year)
{//1 if this is a leap year, else 0
	return (year % 400 == 0 || (year % 4 == 0 && year % 100 != 0));
}

cDate& 
cDate::operator-=(const cPeriod& period)
{//Date - Period = Date
	_year -= period._year;
	int monthDiff = _month - period._month; 
	if (monthDiff <= 0) //If the difference between two months is negative we have to add 12
	{					//and subtract 1 year
		_year -= 1;
		_month = monthDiff + 12;
	}
	else 
		_month = monthDiff;
	int dayDiff = _day - period._day;
	if (dayDiff <= 0) //Case where the difference between days is negative
	{
		if (_month != 1)
			_month -= 1;
		else //If the month is January we have to set December and subtract 1 year
		{
			_month = 12;
			_year -= 1;
		} //Then we have to add a number of days according to the month
		if (_month == 4 || _month == 6 || _month == 9 || _month == 11) //Months with 30 days
			_day = dayDiff + 30;
		else if (_month == 2) //February
			if (IsLeapYear(_year))
				_day = dayDiff + 29;
			else
				_day = dayDiff + 28;
		else //Months with 31 days
			_day = dayDiff + 31;
	}
	else //Case where the difference between days is positive, there is nothing more to do
		_day = dayDiff;
	return *this;
}

cDate 
cDate::operator-(const cPeriod& period) const
{
	cDate res(*this);
	res -= period;
	return res;
}

cDate&
cDate::operator+=(const cPeriod& period)
{//Date + Period = Date
	_year += period._year;
	int monthAdd = _month + period._month;
	if (monthAdd > 12) //If the addition of 2 months is larger than 12 we have to subtract 12
	{				   //and add 1 year
		_year += 1;
		_month = monthAdd - 12;
	}
	else
		_month = monthAdd;
	int dayAdd = _day + period._day; //Here we handle the cases where the addition of days is larger
									 //than a certain number, ACCORDING TO THE MONTH
	if ((_month == 4 || _month == 6 || _month == 9 || _month == 11) && dayAdd > 30) //Months with 30 days
	{
		_month += 1;
		_day = dayAdd - 30;
	}
	else if (_month == 2 && (IsLeapYear(_year) && dayAdd > 29)) //February when the year is leap
	{
		_month += 1;
		_day = dayAdd - 29;
	}
	else if (_month == 2 && (!IsLeapYear(_year) && dayAdd > 28)) //February when the year is not leap
	{
		_month += 1;
		_day = dayAdd - 28;
	}
	else if (dayAdd > 31) //Months with 31 days. December is a special case.
	{
		if (_month != 12)
		{
			_month += 1;
		}
		else
		{
			_year += 1;
			_month = 1;
		}
		_day = dayAdd - 31;
	}
	else
		_day = dayAdd;
	return *this;
}

cDate
cDate::operator+(const cPeriod& period) const
{
	cDate res(*this);
	res += period;
	return res;
}

bool
cDate::operator<(const cDate& d1)
{//Comparison between 2 dates
	if (_year < d1._year) {
		return true;
	}
	if (d1._year < _year) {
		return false;
	}
	else {
		if (_month < d1._month) {
			return true;
		}
		if (d1._month < _month) {
			return false;
		}
		else {
			if (_day < d1._day) {
				return true;
			}
			else {
				return false;
			}
		}
	}
}


