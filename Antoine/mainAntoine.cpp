#include <iostream>z
#include <fstream>
#include <string>
#include <map>

using namespace std;

double conversion(string maturite) {//Conversion of "1Y6M" to 1.5
	int i = 0;
	string s = "";
	double ttm = 0;// time to maturity;
	while (maturite[i] != '\0') {
		while (int(maturite[i]) >= 48 && int(maturite[i]) <= 57) {
			s.push_back(maturite[i]);
			i++;
		}
		if (maturite[i] == 'Y') {
			ttm += stod(s);
			s = "";
			i++;
		}
		if (maturite[i] == 'M') {
			ttm += stod(s) /12;
			s = "";
			i++;
		}
		if (maturite[i] == 'W') {
			ttm += stod(s) * 7 / 360;
			s = "";
			i++;
		}
		if (maturite[i] == 'D') {
			ttm += stod(s) / 360;
			s = "";
			i++;
		}
	}
	return ttm;
}

void replace_coma(string& s) {//Conversion of 1,43 to 1.43
	int i = 0;
	while (s[i] != '\0') {
		if (s[i] == ',') {
			s[i] = '.';
		}
		i++;
	}
}

double interpolation(double T, map<double, double> courbe, string convention = "linear") {//Interpolation for the yield curve from the data
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

map<double, double> courbe_de_taux(string nom) {//Get the market data from a csv - mettre en argument le fichier
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
			replace_coma(taux);
			courbe[conversion(maturite)] = stod(taux)/100;
			maturite = "";
			taux = "";
		}
		data.close();
	}
	else {
		std::cout << "Impossible d'ouvrir le fichier !" << std::endl;
	}
	return courbe;
}

map <double, double> courbe_de_spread(string nom) {
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
			replace_coma(sp);
			spread[conversion(maturite)] = stod(sp)/10000;
			maturite = "";
			sp = "";
		}
		data.close();
	}
	else {
		cout << "Impossible d'ouvrir le fichier !" << endl;
	}
	return spread;
}

double Pricer_IRS(double nominal, double T, map<double, double> courbe, double taux_fixe, double frequence) {//Price an IRS today with same rate for discount and float leg
	double prix=0;
	prix = ZC(interpolation(0,courbe), 0) - ZC(interpolation(T,courbe), T);
	for (double i = frequence; i <= T; i = i + frequence) {
		prix -= ZC(interpolation(i, courbe), i) * frequence * taux_fixe;
	}
	return prix * nominal;
}

double Pricer_IRS_bicourbe(double nominal, double T, map <double, double> discount, map <double, double> ibor, double taux_fixe, double frequence) {
	double prix=0;
	for (double i = frequence; i <= T; i = i + frequence) {
		prix += ZC(interpolation(i, discount), i) * frequence \
			* ((ZC(interpolation(i - frequence, ibor), i - frequence) / ZC(interpolation(i, ibor), i) - 1) / frequence - taux_fixe);
	}
	return prix * nominal;
}

void affiche(map<double, double> mymap) {//show the content of a map
	for (map<double, double>::iterator it = mymap.begin(); it != mymap.end(); it++) {
		std::cout << it->first << " : " << it->second << endl;
	}
}

map <double, double> proba_defaut(map <double, double> spread, map<double, double> discount, double R = 0.4) {
	map <double, double> proba_defaut;
	proba_defaut[0] = 1;
	double L = 1 - R;
	int i = 0;
	for (map<double, double>::iterator it = spread.begin(); it != spread.end(); it++) {
		double proba = 0;
		if (proba_defaut.size() == 1) {
			proba = L / (L + it->first* it->second);
			/*double accru = 40.0 / 360 * it->second;
			cout << accru << endl;
			accru = accru / (ZC(interpolation(it->first, discount), it->first) * (it->first * it->second + L));
			cout << accru << endl;
			proba = proba - accru;*/
		}
		else {
			map<double, double>::iterator it2 = proba_defaut.begin();
			it2++;
			map <double, double>::iterator it4;
			for (it2; it2 != proba_defaut.end(); it2++) {
				it4 = it2;
				it4--;
				proba += ZC(interpolation(it2->first, discount), it2->first) * (L * it4->second - it2->second * (L + (it2->first-it4->first) * it->second));
				it4++;
				if (i == 30) {
					return proba_defaut;
				}
			}
			proba = proba / (ZC(interpolation(it->first, discount), it->first) * (L + (it->first-it4->first) * it->second));
			map<double, double>::iterator it3 = proba_defaut.end();
			it3--;
			proba += L * it3->second / (L + (it->first-it4->first) * it->second);
		}
		i++;
		double accru = 40.0 / 360 * it->second;
		accru = accru / (ZC(interpolation(it->first, discount), it->first) * (it->first * it->second + L));
		proba = proba - accru;
		proba_defaut[it->first] = proba;
	}
	for (map<double, double>::iterator it = proba_defaut.begin(); it != proba_defaut.end(); it++) {
		proba_defaut[it->first] = 1 - it->second;
	}
	return proba_defaut;
}

int main() {
	map<double, double> discount = courbe_de_taux("CourbeOIS06.05.2019");
	map<double, double> libor = courbe_de_taux("CourbeLibor06.05.2019");
	//affiche(discount);
	double taux_fixe = (2.0 / 100);
	//double prix = Pricer_IRS(1, 5, discount, taux_fixe, 0.5);
	//cout << "Le prix avec le mono-courbe est: " << prix << " euros" << endl;
	//double prix2 = Pricer_IRS_bicourbe(1, 5, discount, discount, taux_fixe, 0.5);
	//cout << "Le prix avec le bi-courbe est " << prix2 << " euros" << endl;
	map<double, double> spread = courbe_de_spread("CDSSpreadV");
	affiche(spread);
	map<double, double> prob = proba_defaut(spread, discount);
	//affiche(spread);
	affiche(prob);
	/*map<double, double>::iterator it = spread.begin();
	it++;
	for (it; it != spread.end(); it++) {
		cout << it->first << endl;
	}*/
}