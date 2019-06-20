#define _USE_MATH_DEFINES 
#include "Functions.hh"
#include "cFloatingLeg.hh"

using namespace std;


vector<cDate>
Schedule(const cDate& start,
		 const cPeriod& maturity, 
	     const cPeriod& freq,
		 const eConvention_NonBusinessDay& NonBusinessDayConvention)
{
	if (!start.IsValid() || !maturity.IsValid() ||
		(freq.GetDay() == 0 && freq.GetMonth() == 0 && freq.GetYear() == 0))
	{
		cerr << "Starting date, maturity or frequence is not valid" << endl;
		exit(1);
	}
	cPeriod oneDay(1, 0, 0, conv_ACT_ACT);
	cPeriod twoDay(2, 0, 0, conv_ACT_ACT);
	vector<cDate> schedule;
	cDate tempo = start + maturity;
	if (tempo.WhatDayIsIt() == 6)
	{
		if ( NonBusinessDayConvention == GoForward)
		{
			schedule.push_back(tempo + twoDay);
		}
		if (NonBusinessDayConvention == GoBackward)
		{
			schedule.push_back(tempo - oneDay);
		}
		if (NonBusinessDayConvention == GoToTheClosest)
		{
			schedule.push_back(tempo - oneDay);
		}
	}
	if (tempo.WhatDayIsIt() == 7)
	{
		if (NonBusinessDayConvention == GoForward)
		{
			schedule.push_back(tempo + oneDay);
		}
		if (NonBusinessDayConvention == GoBackward)
		{
			schedule.push_back(tempo - twoDay);
		}
		if (NonBusinessDayConvention == GoToTheClosest)
		{
			schedule.push_back(tempo + oneDay);
		}
	}
	if (tempo.WhatDayIsIt() < 6)
	{
		schedule.push_back(tempo);
	}
	if (maturity.GetDay() == 0 && maturity.GetMonth() == 0 && maturity.GetYear() == 0)
	{
		return schedule;
	}
	while (start < tempo - freq)
	{
		tempo = tempo - freq;
		if (tempo.WhatDayIsIt() == 6)
		{
			if (NonBusinessDayConvention == GoForward)
			{
				schedule.insert(schedule.begin(), tempo + twoDay);
			}
			if (NonBusinessDayConvention == GoBackward)
			{
				schedule.insert(schedule.begin(), tempo-oneDay);
			}
			if (NonBusinessDayConvention == GoToTheClosest)
			{
				schedule.insert(schedule.begin(), tempo-oneDay);
			}
		}
		if (tempo.WhatDayIsIt() == 7)
		{
			if (NonBusinessDayConvention == GoForward)
			{
				schedule.insert(schedule.begin(), tempo+oneDay);
			}
			if (NonBusinessDayConvention == GoBackward)
			{
				schedule.insert(schedule.begin(), tempo-twoDay);
			}
			if (NonBusinessDayConvention == GoToTheClosest)
			{
				schedule.insert(schedule.begin(), tempo+oneDay);
			}
		}
		if (tempo.WhatDayIsIt() < 6)
		{
			schedule.insert(schedule.begin(), tempo);
		}
	}
	return schedule; 
}

void
ShowSchedule(const vector<cDate>& schedule)
{
	for (size_t i = 0; i < schedule.size(); i++)
	{
		schedule[i].Show();
	}
}

cDate 
NumberOfDays_To_Date(int ndays)
{//Convert a number of days since 1/1/1900 to a date
	ndays -= 36526;
	if (ndays <= 0)
	{
		cerr << "Please enter a number of days larger than 36526" << endl;
		exit(1);
	}
	int year = 2000; 
	int month = 1; 
	int day = 1;
	while ((ndays > 366 && cDate::IsLeapYear(year)) || (ndays > 365 && !cDate::IsLeapYear(year)))
	{
		if (ndays > 366 && cDate::IsLeapYear(year))
		{
			ndays -= 366;
			year++;
		}
		else
		{
			ndays -= 365;
			year++;
		}
	}
	vector<int> ndaysPerMonthReverse{ 31,28,31,30,31,30,31,31,30,31,30,31 };
	if (cDate::IsLeapYear(year))
		ndaysPerMonthReverse[1] = 29;
	int i = 0;
	while (ndays > 0)
	{
		if (ndays >= ndaysPerMonthReverse[i])
		{
			ndays -= ndaysPerMonthReverse[i];
			month++;
		}
		else
		{
			day += ndays;
			ndays -= ndaysPerMonthReverse[i];
		}
		i++;
	}
	return cDate(day, month, year);
}

int 
Date_To_NumberOfDays(const cDate& date)
{//Convert a number of days since 1/1/1900 to a date
	int day = date.GetDay();
	int month = date.GetMonth();
	int year = date.GetYear();
	int res = 0;
	for (int countYear = 1900; countYear < year; countYear++)
	{
		if (cDate::IsLeapYear(countYear)) 
			res += 366;
		else 
			res += 365;
	}
	for (int countMonth = 1; countMonth < month; countMonth++)
	{
		if (countMonth == 4 || countMonth == 6 || countMonth == 9 || countMonth == 11) {
			res += 30;
			continue;
		}
		if (countMonth == 2 && cDate::IsLeapYear(year)) {
			res += 29;
			continue;
		}
		if (countMonth == 2) {
			res += 28;
			continue;
		}
		else
			res += 31;
	}
	res += day;
	return res + 1;
}
