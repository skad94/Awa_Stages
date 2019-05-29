#include "cDate.hh"
#include "cPeriod.hh"
#include "Functions.hh"


int main()
{
	/*cDate d1(1, 1, 1900);
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
	d1.Show();*/
	cDate start(18, 6, 2019);
	cPeriod maturity(0, 0, 5, { conv_30_360 });
	cPeriod freq(0, 6, 0, { conv_30_360 });
	vector <cDate> schedule = Schedule(start,maturity,freq);
	ShowSchedule(schedule);
	cout << endl;
	NumberOfDays_To_Date(45760).Show();
	//cout << cDate::IsLeapYear(1900);
}