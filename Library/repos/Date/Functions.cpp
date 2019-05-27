#define _USE_MATH_DEFINES 
#include "Functions.hh"

using namespace std;


vector<cDate>
Schedule(cDate& start,const cPeriod& maturity, const cPeriod& freq) 
{
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
