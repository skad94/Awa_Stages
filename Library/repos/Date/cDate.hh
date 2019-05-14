
class cDate 
{
private:
	int _day;
	int _month;
	int _year;
public:
	cDate();
	cDate(const int& day, const int& month, const int& year);
	cDate(const cDate& d);
	~cDate();
	void ShowDate();
};



