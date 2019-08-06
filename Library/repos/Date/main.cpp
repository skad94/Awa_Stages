#include "cDate.hh"
#include "cPeriod.hh"
#include "cFloatingLeg.hh"
#include "cFixedLeg.hh"
#include "cInterestRateSwap.hh"
#include <iostream>
using namespace std;

void primes(int N)
{
	for (int i = 2; i <= N; i++)
	{
		int j;
		for (j = 2; j < (int) sqrt(i) + 1; j++)
		{
			if (i % j == 0) break;
		}
		if (j == (int) sqrt(i) + 1) cout << i << endl;
	}
}



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

	eConvention conv = conv_30_360;
	cDate start(28, 6, 2019);
	cPeriod maturity(0, 0, 5, conv);
	cPeriod freq(0, 6, 0, conv);
	//vector <cDate> schedule = Schedule(start, maturity, freq, GoForward);
	//ShowSchedule(schedule);
	//cout << maturity.ConvertToDayFraction() << endl;
	//cout << freq.ConvertToDayFraction() << endl;
	//(cDate(1, 1, 2018).minus(cDate(2, 3, 2015), conv_30_360)).Show();
	//(cDate(29, 8, 2019) + cPeriod(0, 6, 0, conv_30_360)).Show();*/
	cPeriod paymentGap(0, 0, 0, conv);
	cPeriod fixingGap(0, 0, 0, conv);

	map<double, double> discount = YieldCurve("CourbeOIS28.06.2019", conv, "lroussel");
	double fixedRate = 3.0 / 100;
	cLeg leg(GoForward, conv, maturity, freq, start, paymentGap, discount);
	const cFixedLeg fixedLeg(leg, fixedRate);
	const cFloatingLeg floatingLeg(leg, fixingGap);
	cInterestRateSwap swap(fixedLeg, floatingLeg);
	double notional = 70000;
	cout << floatingLeg.PriceLeg() * notional << endl;
	cout << fixedLeg.PriceLeg() * notional << endl;
	cout << swap.Price_IRS() * notional << endl;
	primes(100);
	//double frac = (cDate(3, 7, 2023).minus(cDate(27, 6, 2019), conv_30_360)).ConvertToDayFraction();
	//cout << ZC(Interpolation(frac, discount), frac) << endl;
}