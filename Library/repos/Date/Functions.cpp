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

/*cDate 
NumberOfDays_To_Date(int ndays)
{
	int year = 1900; 
	int month = 1; 
	int day = 1;
	while (ndays > 0) 
	{
		if (ndays > 366 && cDate::IsLeapYear(year))
		{
			ndays -= 366;
			year++;
		}
		else if (ndays > 365 && !cDate::IsLeapYear(year))
		{
			
		}
			if (ndays > 366)
			{
				year++;
				ndays -= 366;
			}
			else
			{


		}
		else
		{
			annee++;
			ndays - 365
		}

}*/
