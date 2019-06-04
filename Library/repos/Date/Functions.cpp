#define _USE_MATH_DEFINES 
#include "Functions.hh"

using namespace std;


vector<cDate>
Schedule(cDate& start, const cPeriod& maturity, const cPeriod& freq)
{
	if (!start.IsValid() || !maturity.IsValid() || (freq.GetDay() == 0 && freq.GetMonth() == 0 && freq.GetYear() == 0))
	{
		cerr << "Starting date, maturity or frequence is not valid" << endl;
		exit(1);
	}
	vector<cDate> schedule;
	cDate tempo = start + maturity;
	if (maturity.GetDay() == 0 && maturity.GetMonth() == 0 && maturity.GetYear() == 0)
	{
		schedule.push_back(tempo);
		return schedule;
	}
	while (start < tempo - freq)
	{
		tempo = tempo - freq;
		schedule.insert(schedule.begin(), tempo);
	}
	return schedule;
}

void
ShowSchedule(vector<cDate>& schedule)
{
	int size = schedule.size();
	for (int i = 0; i < size ; i++)
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
		cerr << "Please enter a date after 1/1/2000" << endl;
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
