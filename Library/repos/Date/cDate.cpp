#include "cDate.hh"
#include "cPeriod.hh"

cDate::cDate() 
	:_day(0), _month(0), _year(0) {}

cDate::cDate(const int& day, const int& month, const int& year)
	: _day(day), _month(month), _year(year) 
{
	if (!IsValid()) 
	{
		SetAsInvalid();
	}
}

cDate::cDate(const cDate& d) 
	: _day(d._day), _month(d._month), _year(d._year) 
{
	if (!IsValid())
	{
		SetAsInvalid();
	}
}

void 
cDate::ShowDate()
{
	cout << _day << "/" << _month << "/" << _year << endl;
}

bool 
cDate::IsValid() const
{//Checking whether the date is valid. Values can be 0 since the date can be a period.
	if (0 <= _month && _month <= 12)
	{
		if (_month == 4 || _month == 6 || _month == 9 || _month == 11) //Months with 30 days
			return (0 <= _day && _day <= 30);
		else if (_month == 2) //February
			if (IsLeapYear(_year))
				return (0 <= _day && _day <= 29);
			else
				return (0 <= _day && _day <= 28);
		else //Months with 31 days or a period with _month = 0
			return (0 <= _day && _day <= 31);
	}
	return 0;
}

void 
cDate::SetAsInvalid()
{
	_day = 0;
	_month = 0;
	_year = 0;
}

bool 
cDate::IsLeapYear(const int& year)
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
			if (IsLeapYear(_year))
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
cDate::operator-(const cPeriod& period) const
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
	else if (_month == 2 && (IsLeapYear(_year) && dayAdd > 29))
	{
		_month += 1;
		_day = dayAdd - 29;
	}
	else if (_month == 2 && (!IsLeapYear(_year) && dayAdd > 28))
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
cDate::operator+(const cPeriod& period) const
{
	cDate res(*this);
	res += period;
	return res;
}

/*bool 
cDate::operator<(const Date& d1) const 
{
	if (annee < d1.annee) {
		return true;
	}
	if (d1.annee < annee) {
		return false;
	}
	else {
		if (mois < d1.mois) {
			return true;
		}
		if (d1.mois < mois) {
			return false;
		}
		else {
			if (jour < d1.jour) {
				return true;
			}
			else {
				return false;
			}
		}
	}
}*/

