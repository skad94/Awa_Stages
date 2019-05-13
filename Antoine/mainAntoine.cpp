#include <iostream>z
#include <fstream>
#include <string>
#include <map>

using namespace std;

double ConversionReuters(string maturite) {//Conversion of "1Y6M" to 1.5
	int i = 0;
	string s = "";
	double timetomaturity = 0;
	while (maturite[i] != '\0') {
		while (int(maturite[i]) >= 48 && int(maturite[i]) <= 57) {
			s.push_back(maturite[i]);
			i++;
		}
		if (maturite[i] == 'Y') {
			timetomaturity += stod(s);
			s = "";
			i++;
		}
		if (maturite[i] == 'M') {
			timetomaturity += stod(s) /12;
			s = "";
			i++;
		}
		if (maturite[i] == 'W') {
			timetomaturity += stod(s) * 7 / 360;
			s = "";
			i++;
		}
		if (maturite[i] == 'D') {
			timetomaturity += stod(s) / 360;
			s = "";
			i++;
		}
	}
	return timetomaturity;
}

void ReplaceComa(string& s) {//Conversion of 1,43 to 1.43
	int i = 0;
	while (s[i] != '\0') {
		if (s[i] == ',') {
			s[i] = '.';
		}
		i++;
	}
}

double Interpolation(double T, map<double, double> courbe, string convention = "linear") {//Interpolation for the yield curve from the data
	double t1;
	double t2;
	double r1;
	double r2;
	map<double, double>::iterator it = courbe.begin();
	if (T < it->first) {
		return it->second;
	}
	for (it; it != courbe.end(); it++) {
		t2 = it->first; r2 = it->second;
		if (T < it->first) {
			it--;
			t1 = it->first; r1 = it->second;
			return (T - t2) / (t1 - t2) * r1 + (T - t1) / (t2 - t1) * r2;
		}
	}
	return r2;
}

double ZC(double taux, double T, string convention = "compose", double t = 0) {//ZC Price with 2 convention T the maturity and t the date of valuation are exprimed in years
	if (convention == "compose")
		return pow(1 / (1 + taux), T - t);
		//return 1;
	if (convention == "expo")
		return exp(-taux * (T - t));
}

map<double, double> YieldCurve(string nom) {//Get the market data from a csv - mettre en argument le fichier
	map<double, double> courbe;
	std::ifstream data("C:/users/alepeltier/Documents/Data/"+nom+".csv", std::ios::in);
	if (data) {
		string contenu;
		string maturite= "";
		string taux = "";
		for (int i = 0; i < 6; i++) {
			getline(data, contenu);
		}
		while (getline(data, contenu)) {
			int i = 0;
			while (contenu[i] != ';') {
				maturite.push_back(contenu[i]);
				i++;
			}
			int nb = 0;
			while (nb < 3) {
				if (contenu[i] != ';') {
					i++;
				}
				else {
					nb++; i++;
				}
			}
			while (contenu[i] != ';') {
				taux.push_back(contenu[i]);
				i++;
			}
			ReplaceComa(taux);
			courbe[ConversionReuters(maturite)] = stod(taux)/100;
			maturite = "";
			taux = "";
		}
		data.close();
	}
	else {
		std::cout << "Can't open the file !" << std::endl;
	}
	return courbe;
}

map <double, double> SpreadCurve(string nom) {
	map <double, double> spread;
	std::ifstream data("C:/users/alepeltier/Documents/Data/" + nom + ".csv", std::ios::in);
	if (data) {
		string contenu;
		string maturite = "";
		string sp = "";
		for (int i = 0; i < 1; i++) {
			getline(data, contenu);

		}
		while (getline(data, contenu)) {
			int i = 0;
			while (contenu[i] != ';') {
				maturite.push_back(contenu[i]);
				i++;
			}
			int nb = 0;
			while (nb < 2) {
				if (contenu[i] != ';') {
					i++;
				}
				else {
					nb++; i++;
				}
			}
			while (contenu[i] != ';') {
				sp.push_back(contenu[i]);
				i++;
			}
			cout << sp << endl;
			ReplaceComa(sp);
			spread[ConversionReuters(maturite)] = stod(sp)/10000;
			maturite = "";
			sp = "";
		}
		data.close();
	}
	else {
		cout << "Can't open the file !" << endl;
	}
	return spread;
}

double PricerIRS(double nominal, double T, map<double, double> discount, double fix_rate, double frequency) {//Price an IRS today with same rate for discount and float leg
	double prix=0;
	prix = ZC(Interpolation(0,discount), 0) - ZC(Interpolation(T,discount), T);
	for (double i = frequency; i <= T; i = i + frequency) {
		prix -= ZC(Interpolation(i, discount), i) * frequency * fix_rate;
	}
	return prix * nominal;
}

double PricerIRSBicourbe(double nominal, double T, map <double, double> discount, map <double, double> ibor, double fix_rate, double frequency) {
	double prix=0;
	for (double i = frequency; i <= T; i = i + frequency) {
		prix += ZC(Interpolation(i, discount), i) * frequency \
			* ((ZC(Interpolation(i - frequency, ibor), i - frequency) / ZC(Interpolation(i, ibor), i) - 1) / frequency - fix_rate);
	}
	return prix * nominal;
}

void Affiche(map<double, double> mymap) {//show the content of a map
	for (map<double, double>::iterator it = mymap.begin(); it != mymap.end(); it++) {
		std::cout << it->first << " : " << it->second << endl;
	}
}

map <double, double> DefaultProb(map <double, double> spread, map<double, double> discount, double R = 0.4) {
	map <double, double> proba_defaut;
	proba_defaut[0] = 1;
	double L = 1 - R;//LOSS GIVEN
	int i = 0;
	for (map<double, double>::iterator it = spread.begin(); it != spread.end(); it++) {
		double proba = 0;
		if (proba_defaut.size() == 1) {
			proba = L / (L + it->first* it->second);
		}
		else {
			map<double, double>::iterator it2 = proba_defaut.begin();
			it2++;
			map <double, double>::iterator it4;
			for (it2; it2 != proba_defaut.end(); it2++) {
				it4 = it2;
				it4--;
				proba += ZC(Interpolation(it2->first, discount), it2->first) * (L * it4->second - it2->second * (L + (it2->first-it4->first) * it->second));
				it4++;
				if (i == 30) {
					return proba_defaut;
				}
			}
			proba = proba / (ZC(Interpolation(it->first, discount), it->first) * (L + (it->first-it4->first) * it->second));
			map<double, double>::iterator it3 = proba_defaut.end();
			it3--;
			proba += L * it3->second / (L + (it->first-it4->first) * it->second);
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

int main() {
	map<double, double> discount = YieldCurve("CourbeOIS06.05.2019");
	map<double, double> libor = YieldCurve("CourbeLibor06.05.2019");
	//Affiche(discount);
	double fix_rate = (2.0 / 100);
	//double prix = PricerIRS(1, 5, discount, fix_rate, 0.5);
	//cout << "Le prix avec le mono-courbe est: " << prix << " euros" << endl;
	//double prix2 = PricerIRSBicourbe(1, 5, discount, discount, fix_rate, 0.5);
	//cout << "Le prix avec le bi-courbe est " << prix2 << " euros" << endl;
	map<double, double> spread = SpreadCurve("CDSSpreadV");
	Affiche(spread);
	map<double, double> prob = DefaultProb(spread, discount);
	//Affiche(spread);
	Affiche(prob);
	/*map<double, double>::iterator it = spread.begin();
	it++;
	for (it; it != spread.end(); it++) {
		cout << it->first << endl;
	}*/
}