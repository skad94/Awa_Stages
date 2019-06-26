#include "cDate.hh"
#include "cPeriod.hh"
#include "cFloatingLeg.hh"
#include "cFixedLeg.hh"
#include <iostream>
using namespace std;



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
	cDate start(1, 9, 2019);
	cPeriod maturity(0, 0, 1, conv_30_360);
	cPeriod freq(0, 6, 0, conv_30_360);
	vector <cDate> schedule = Schedule(start, maturity, freq, GoForward);
	ShowSchedule(schedule);
	cout << maturity.ConvertToDayFraction() << endl;
	cout << freq.ConvertToDayFraction() << endl;
	(cDate(1, 1, 2018).minus(cDate(2, 3, 2015), conv_30_360)).Show();
}