#include "cDate.hh"
#include "cPeriod.hh"

vector<cDate>
Schedule(cDate& date, const cPeriod& maturity, const cPeriod& freq)
{
	vector <cDate> schedule;
	cDate tempo = date + maturity;
	schedule.push_back(tempo);
	while 
}

int main() 
{
	cDate d(1, 2, 1959);
	cPeriod p(30, 0, 9, 30, 360);
	(d -= p).ShowDate();
	cDate tempo = d+p;
	tempo.ShowDate();

	return 0;
}