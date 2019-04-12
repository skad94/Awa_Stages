import numpy as np
import os
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

os.chdir('C:\\Users\lroussel\Documents\Données Reuters') # à modifier

df = pd.read_excel('CAC 40 vol daily.xlsx', index_col=None, header=None) # lit le fichier excel
aa = np.flip(np.array(df)[1:])
n = len(aa)

def m(q,delta): # fonction qui renvoie le moment d'ordre q des différences de log-volatilité en fonction du lag delta
    diff_log = []
    i = 0
    while i + delta < n:
        diff_log.append(np.log(float(aa[i + delta]/aa[i])))
        i += delta
    res = np.abs(np.array(diff_log))**q
    res = np.mean(res)
    return res

# On trace log(m(q,delta)) en fonction de log(delta), pour différentes valeurs de q
N = 15 # nombre de valeurs de q
q_ = np.arange(1,N+1)/5
x = np.log(np.arange(1,int(np.exp(3)))) # x = log(delta)
y = [[np.log(m(q,delta)) for delta in np.arange(1,int(np.exp(3)))] for q in q_]
for i in y: 
    plt.plot(x,i)
plt.xlabel("log(delta)")
plt.ylabel("log(m(q,delta))")
plt.title("S&P 500 : log(m(q,delta)) en fonction de log(delta)")
plt.legend(["q = "+ str(q) for q in q_])
plt.show()

# On fait une régression linéaire pour chaque valeur de q
x_ = x.reshape(-1,1)
reg = [LinearRegression().fit(x_,y[i]) for i in range(N)]
scores = [reg[i].score(x_,y[i]) for i in range(N)]
print("scores = " + str(scores))
coefs = [reg[i].coef_ for i in range(N)] # on obtient la pente pour chaque q
intercept = [reg[i].intercept_ for i in range(N)]
print("intercept = " + str(intercept))

# Régression linéaire sur la pente en fonction de q
# Le coefficient directeur qu'on va obtenir correspond à H le paramètre de Hurst
reg_pente = LinearRegression().fit(q_.reshape(-1,1),coefs)
plt.plot(q_,coefs,'o')
plt.plot(q_,np.dot(q_.reshape(-1,1),reg_pente.coef_) + reg_pente.intercept_)
plt.xlabel("q")
plt.ylabel("coef")
plt.title("S&P 500 : pente en fonction de q. On obtient H = " + str(float(reg_pente.coef_)))
plt.legend(['coef','linear fit'])
plt.show()

