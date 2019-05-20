#include "cDate.hh"
#include "cPeriod.hh"

int main() 
{
	cDate d(1, 2, 1959);
	cPeriod p(30, 0, 9, 30, 360);
	(d -= p).ShowDate();



	return 0;
}