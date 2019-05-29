#define _USE_MATH_DEFINES 
#include "Functions.hh"

using namespace std;


vector<cDate>
Schedule(cDate& start,const cPeriod& maturity, const cPeriod& freq) 
{
	if (!start.IsValid() || !maturity.IsValid())
	{
		cerr << "Starting date or maturity is not valid" << endl;
		exit(1);
	}
	vector<cDate> schedule;
	cDate tempo = start + maturity ;
	schedule.push_back(tempo);
	while (start < tempo - freq )
	{
		tempo = tempo - freq;
		schedule.insert( schedule.begin() ,tempo);
	}
	schedule.insert(schedule.begin(), start);
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
	int year = 1900; 
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
	vector<int> ndaysPerMonthReverse{ 31,30,31,30,31,31,30,31,30,31,28,31 };
	if (cDate::IsLeapYear(year))
		ndaysPerMonthReverse[10] = 29;
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
	}
}
