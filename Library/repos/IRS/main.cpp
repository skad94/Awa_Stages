#include <iostream>
#include <fstream>
#include <string>
#include <map>

using namespace std;

double
ConversionReuters(const string maturity)
{//Conversion of "1Y6M" to 1.5
	int i = 0;
	string s = "";
	double timeToMaturity = 0;
	while (maturity[i] != '\0') {
		while (int(maturity[i]) >= 48 && int(maturity[i]) <= 57)
		{
			s.push_back(maturity[i]);
			i++;
		}
		if (maturity[i] == 'Y')
		{
			timeToMaturity += stod(s);
			s = "";
			i++;
		}
		if (maturity[i] == 'M')
		{
			timeToMaturity += stod(s) / 12;
			s = "";
			i++;
		}
		if (maturity[i] == 'W')
		{
			timeToMaturity += stod(s) * 7 / 360;
			s = "";
			i++;
		}
		if (maturity[i] == 'D')
		{
			timeToMaturity += stod(s) / 360;
			s = "";
			i++;
		}
	}
	return timeToMaturity;
}

void
ReplaceComa(string& s)
{//Conversion of 1,43 to 1.43
	int i = 0;
	while (s[i] != '\0')
	{
		if (s[i] == ',')
		{
			s[i] = '.';
		}
		i++;
	}
}

double
Interpolation(const double& x,
	const map<double, double> curve,
	const string convention = "linear")
{//Interpolation for the yield curve from the data
	double t1;
	double t2;
	double r1;
	double r2;
	map<double, double>::const_iterator it = curve.begin();
	if (x < it->first)
	{
		return it->second;
	}
	for (it; it != curve.end(); it++)
	{
		t2 = it->first; r2 = it->second;
		if (x < it->first)
		{
			it--;
			t1 = it->first; r1 = it->second;
			return (x - t2) / (t1 - t2) * r1 + (x - t1) / (t2 - t1) * r2;
		}
	}
	return r2;
}

double
ZC(const double& rate,
	const double& maturity,
	string convention = "composed",
	const double& t = 0)
{//ZC Price with 2 convention T the maturity and t the date of valuation are exprimed in years
	if (convention == "composed")
	{
		return pow(1.0 / (1.0 + rate), maturity - t);
		//return 1;
	}
	if (convention == "expo")
	{
		return exp(-rate * (maturity - t));
	}
	else {
		return 0;
	}
}

map<double, double>
YieldCurve(string nom,
	string user = "alepeltier")
{//Get the market data from a csv - mettre en argument le fichier
	map<double, double> courbe;
	std::ifstream data("C:/users/" + user + "/Documents/GitHub/Awa_Stages/Data/" + nom + ".csv", std::ios::in);
	if (data)
	{
		string contenu;
		string maturite = "";
		string taux = "";
		for (int i = 0; i < 6; i++)
		{
			getline(data, contenu);
		}
		while (getline(data, contenu))
		{
			int i = 0;
			while (contenu[i] != ';')
			{
				maturite.push_back(contenu[i]);
				i++;
			}
			int nb = 0;
			while (nb < 3)
			{
				if (contenu[i] != ';')
				{
					i++;
				}
				else
				{
					nb++; i++;
				}
			}
			while (contenu[i] != ';')
			{
				taux.push_back(contenu[i]);
				i++;
			}
			ReplaceComa(taux);
			courbe[ConversionReuters(maturite)] = stod(taux) / 100;
			maturite = "";
			taux = "";
		}
		data.close();
	}
	else
	{
		std::cout << "Can't open the file !" << std::endl;
	}
	return courbe;
}

map <double, double>
SpreadCurve(string nom,
	string user = "alepeltier")
{
	map <double, double> spread;
	std::ifstream data("C:/users/" + user + "/Documents/GitHub/Awa_Stages/Data/" + nom + ".csv", std::ios::in);
	if (data)
	{
		string contenu;
		string maturite = "";
		string sp = "";
		for (int i = 0; i < 1; i++)
		{
			getline(data, contenu);

		}
		while (getline(data, contenu))
		{
			int i = 0;
			while (contenu[i] != ';')
			{
				maturite.push_back(contenu[i]);
				i++;
			}
			int nb = 0;
			while (nb < 2)
			{
				if (contenu[i] != ';')
				{
					i++;
				}
				else
				{
					nb++; i++;
				}
			}
			while (contenu[i] != ';')
			{
				sp.push_back(contenu[i]);
				i++;
			}
			ReplaceComa(sp);
			spread[ConversionReuters(maturite)] = stod(sp) / 10000;
			maturite = "";
			sp = "";
		}
		data.close();
	}
	else
	{
		cout << "Can't open the file !" << endl;
	}
	return spread;
}

double
PricerIRS(const double& nominal,
	const double maturity,
	map<double, double> discount,
	double fix_rate,
	double frequency)
{//Price an IRS today with same rate for discount and float leg
	double prix = 0;
	prix = ZC(Interpolation(0, discount), 0) - ZC(Interpolation(maturity, discount), maturity);
	for (double i = frequency; i <= maturity; i = i + frequency)
	{
		prix -= ZC(Interpolation(i, discount), i) * frequency * fix_rate;
	}
	return prix * nominal;
}

double
PricerIRSBicourbe(const double& nominal,
	const double& maturity,
	const map <double, double> discount,
	const map <double, double> ibor,
	double& fix_rate,
	double& frequency)
{
	double prix = 0;
	for (double i = frequency; i <= maturity; i = i + frequency)
	{
		double liborForwardRate = (1 / frequency) \
			* (ZC(Interpolation(i - frequency, ibor), i - frequency) / ZC(Interpolation(i, ibor), i) - 1);
		prix += ZC(Interpolation(i, discount), i) * frequency \
			* (liborForwardRate - fix_rate);
	}
	return prix * nominal;
}

void
Affiche(const map<double, double> mymap)
{//show the content of a map
	for (map<double, double>::const_iterator it = mymap.begin(); it != mymap.end(); it++)
	{
		std::cout << it->first << " : " << it->second << endl;
	}
}

map <double, double>
DefaultProb(const map <double, double> spread,
	const map<double, double> discount,
	const double& recoveryRate = 0.4)
{
	map <double, double> proba_defaut;
	proba_defaut[0] = 1;
	double L = 1 - recoveryRate;//LOSS GIVEN
	int i = 0;
	for (map<double, double>::const_iterator it = spread.begin(); it != spread.end(); it++)
	{
		double proba = 0;
		if (proba_defaut.size() == 1)
		{
			proba = L / (L + it->first * it->second);
		}
		else
		{
			map<double, double>::const_iterator it2 = proba_defaut.begin();
			it2++;
			map <double, double>::const_iterator it4;
			for (it2; it2 != proba_defaut.end(); it2++)
			{
				it4 = it2;
				it4--;
				proba += ZC(Interpolation(it2->first, discount), it2->first) \
					* (L * it4->second - it2->second * (L + (it2->first - it4->first) * it->second));
				it4++;
				if (i == 30)
				{
					return proba_defaut;
				}
			}
			proba = proba / (ZC(Interpolation(it->first, discount), it->first) \
				* (L + (it->first - it4->first) * it->second));
			map<double, double>::const_iterator it3 = proba_defaut.end();
			it3--;
			proba += L * it3->second / (L + (it->first - it4->first) * it->second);
		}
		i++;
		double accru = 40.0 / 360 * it->second;
		accru = accru / (ZC(Interpolation(it->first, discount), it->first) * (it->first * it->second + L));
		proba = proba - accru;
		proba_defaut[it->first] = proba;
	}
	for (map<double, double>::iterator it = proba_defaut.begin(); it != proba_defaut.end(); it++) {
		proba_defaut[it->first] = 1 - it->second;
	}
	return proba_defaut;
}

void
Output(map<double, double> discount, string user = "alepeltier")
{
	map<double, double> courbe;
	std::ofstream data("C:/users/" + user + "/Documents/GitHub/Awa_Stages/Data/Output.csv", std::ios::trunc);
	if (data)
	{
		for (map<double, double>::const_iterator it = discount.begin(); it != discount.end(); it++)
		{
			data << it->first << "," << it->second << endl;
		}
		data.close();
	}
	else
	{
		std::cout << "Can't open the file !" << std::endl;
	}
}

int main() {
	map<double, double> discount = YieldCurve("CourbeOIS06.05.2019");
	//map<double, double> libor = YieldCurve("CourbeLibor06.05.2019");
	double fix_rate = (2.0 / 100);
	map<double, double> spread = SpreadCurve("CDSSpreadV");
	map<double, double> prob = DefaultProb(spread, discount);
	Affiche(prob);
}