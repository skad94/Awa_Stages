#include "../Date/cDate.hh"
#include "../Date/cPeriod.hh"
#include "../Math/cSquareMatrix.hh"
#include "../Math/Functions.hh"

int main() 
{
	cout << "This is the Test project." << endl;
	cout << "Running tests for cDate and cPeriod ... ";
	cDate d1(1, 1, 1900);
	cPeriod p1(1, 3, 119, conv_30_360);
	d1 += p1;
	cPeriod p2(30, 2, 0, conv_30_360);
	d1 -= p2;
	d1 -= p1;
	d1 += p2;
	cDate d2(18, 0, 2000);
	cPeriod p3(18, 0, 2000, conv_ACT_ACT);
	int d2_valid = d2.IsValid();
	int p3_valid = p3.IsValid();
	if (d1 < d2)
		d2.SetAsInvalid();
	cDate d3(-1, -1, -1);
	int comp = ((d2 < d3) || (d3 < d2));
	int d1_day = d1.GetDay();
	int d1_month = d1.GetMonth();
	int d1_year = d1.GetYear();
	if (d1_day == 1 && d1_month == 1 && d1_year == 1900 && d2_valid == 0 && p3_valid == 1 && comp == 0)
		cout << "PASSED" << endl;
	else
		cout << "FAILED" << endl;
}