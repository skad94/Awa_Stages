import numpy as np
import os
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

os.chdir('C:\\Users\lroussel\Documents\Données Reuters') # à modifier


# Extraction des données :

data_ = 'FOAT_vol'
df = pd.read_excel(data_ + ".xlsx", index_col=None, header=None) # lit le fichier excel
aa = np.flip(np.array(df)[1:])
n = len(aa)


# Etude de l'aspect "rough" grâce à la volatilité réalisée et recherche de l'exposant de Hurst H :

def m(q,delta): # fonction qui renvoie le moment d'ordre q des différences de log-volatilité en fonction du lag delta
    diff_log = []
    i = 0
    while i + delta < n:
        diff_log.append(np.log(float(aa[i + delta]/aa[i])))
        i += delta
    res = np.abs(np.array(diff_log))**q
    res = np.mean(res)
    return res

def regr():
    # On trace log(m(q,delta)) en fonction de log(delta), pour différentes valeurs de q
    N = 4 # nombre de valeurs de q
    q_ = np.arange(1,N+1)/2
    x = np.log(np.arange(1,int(np.exp(3)))) # x = log(delta)
    y = [[np.log(m(q,delta)) for delta in np.arange(1,int(np.exp(3)))] for q in q_]
    plt.figure(1)
    for i in y: 
        plt.plot(x,i,'o')
    plt.xlabel("log(delta)")
    plt.ylabel("log(m(q,delta))")
    plt.title(data_ + " : log(m(q,delta)) en fonction de log(delta)")
    plt.legend(["q = "+ str(q) for q in q_])

    # On fait une régression linéaire pour chaque valeur de q
    x_ = x.reshape(-1,1)
    reg = [LinearRegression().fit(x_,y[i]) for i in range(N)]
    scores = [reg[i].score(x_,y[i]) for i in range(N)]
    print("scores = " + str(scores))
    coefs = [reg[i].coef_ for i in range(N)] # on obtient la pente pour chaque q
    #print(np.exp(reg[N-1].intercept_/2))   # estimation du paramètre "nu"
    
    # Régression linéaire sur la pente en fonction de q
    # Le coefficient directeur qu'on va obtenir correspond à H le paramètre de Hurst
    reg_pente = LinearRegression().fit(q_.reshape(-1,1),coefs)
    plt.figure(2)
    plt.plot(q_,coefs,'o')
    plt.plot(q_,np.dot(q_.reshape(-1,1),reg_pente.coef_) + reg_pente.intercept_)
    plt.xlabel("q")
    plt.ylabel("coef")
    plt.title(data_ + " : pente en fonction de q. On obtient H = " + str(round(float(reg_pente.coef_),2)))
    plt.legend(['coef','linear fit'])
    plt.show()

#regr()



# Ce morceau de code sert à retrouver H à partir d'une simulation trajectoire de MBf :

data_ = 'export_fBMs' 
df = pd.read_excel(data_ + ".xlsx")
aa = np.array(df)[1:]
jj = [aa[:,i] for i in range(9)]
n = len(jj[0])
for i in range(9):
   jj[i] = jj[i] - np.min(jj[i])*2
aa = jj[4] # modifier cette ligne en mettant un entier entre 0 et 8 entre les crochets
           # 0 pour H = 0.1; 1 pour H = 0.2; ... ; 8 pour H = 0.9
#regr()    



# Etude du skew ATM pour SPX :

def regr_skew():
    data_ATM = 'SPX ATM 15042019'
    rr = pd.read_excel(data_ATM + ".xlsx") # lit le fichier excel
    rr = np.array(rr)
    df2 = rr[:,0]
    dates = rr[:,1]
    a = df2[1:]
    a = [float(z) for z in a]
    dates = dates[1:]
    dates = [float(z) for z in dates]
    n_ = len(a)
    skew = []
    for i in range(n_-1): # approximation du skew ATM par différences finies décentrées à droite
        skew.append((a[i+1] - a[i])/(dates[i+1] - dates[i]))

    # Régression linéaire sur log(skew) en fonction de log(maturité) pour trouver l'exposant alpha
    dates = np.log(dates[:n_-1])
    skew = np.log(skew)
    #dates = np.concatenate((dates[:2],dates[3:]))[:10] # FTSE 100
    #skew = np.concatenate((skew[:2],skew[3:]))[:10]    # FTSE 100
    reg_skew = LinearRegression().fit(dates.reshape(-1,1),skew)

    plt.figure(3)
    plt.plot(np.exp(dates),np.exp(skew),'o')
    plt.plot(np.exp(dates),np.exp(np.dot(dates.reshape(-1,1),reg_skew.coef_) + reg_skew.intercept_))
    plt.title("SPX : ATM skew en fonction de la maturité. alpha = " + str(round(-float(reg_skew.coef_),2)))
    plt.legend(['ATM skew','fit'])
    #print(reg_skew.score(dates.reshape(-1,1),skew))
    plt.show()

#regr_skew()