#include "cDate.hh"
#include "cPeriod.hh"

/*vector<cDate>
Schedule(cDate& date, const cPeriod& maturity, const cPeriod& freq)
{
	vector <cDate> schedule;
	cDate tempo = date + maturity;
	schedule.push_back(tempo);
	while 
}*/

int main()
{
	cDate d1(1, 1, 1900);
	d1.Show();
	cPeriod p1(1, 3, 119, conv_30_360);
	p1.Show();
	d1 += p1;
	d1.Show();
	cPeriod p2(30, 2, 0, conv_30_360);
	p2.Show();
	d1 -= p2;
	d1.Show();
	(d1 - p1).Show();
	d1 -= p1;
	d1 += p2;
	d1.Show();
}