
class cDate 
{
private:
	int _day;
	int _month;
	int _year;
public:
	cDate();
	cDate(int day, int month, int year);
	cDate(const cDate& d);
	~cDate();
	void ShowDate();
};



