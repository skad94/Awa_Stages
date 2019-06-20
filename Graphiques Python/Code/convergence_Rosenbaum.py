import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
import time

start_date = time.time()

def simul_Prix(T,beta,phi1,phi2,mu): # On simule le prix microscopique à l'aide de processus de Hawkes
    tau1 = []                        # Les temps de sauts de N+
    tmp1 = 0
    n1 = 0
    while tmp1 < T:
        lambda_barre = mu + np.sum([phi1(tmp1 - t) + beta*phi2(tmp1 - t) for t in tau1])
        u = np.random.uniform()
        w = (-1/lambda_barre)*np.log(u)
        tmp1 += w
        d = np.random.uniform()
        if lambda_barre*d < mu + np.sum([phi1(tmp1 - t) + beta*phi2(tmp1 - t) for t in tau1]):
            n1 += 1
            tau1.append(tmp1)
        print(tmp1)
    if tmp1 > T and tmp1 in tau1:
        tau1.remove(tmp1)
    tau2 = []                        # Les temps de sauts de N-
    tmp2 = 0
    n2 = 0
    while tmp2 < T:
        lambda_barre = mu + np.sum([phi1(tmp2 - t) + beta*phi2(tmp2 - t) for t in tau2])
        u = np.random.uniform()
        w = (-1/lambda_barre)*np.log(u)
        tmp2 += w
        d = np.random.uniform()
        if lambda_barre*d < mu + np.sum([phi1(tmp2 - t) + beta*phi2(tmp2 - t) for t in tau2]):
            n2 += 1
            tau2.append(tmp2)
    if tmp2 > T and tmp2 in tau2:
        tau2.remove(tmp2)
    tau_ = np.unique(np.concatenate((tau1,tau2)))
    n = len(tau_)
    Np = np.zeros(n)                 # N+
    Nm = np.zeros(n)                 # N-
    k = 1                            # Ici on cherche à savoir, pour chaque instant de saut, lequel des deux N+ ou N- saute afin de pouvoir
    i = 1                            # représenter la trajectoire du prix P = N+ - N-
    j = 1
    n1 = len(tau1)
    n2 = len(tau2)
    while i < n1 and j < n2:
        if tau1[i] < tau2[j]:
            Np[k] = Np[k-1] + 1
            Nm[k] = Nm[k-1]
            i += 1
        else:
            Nm[k] = Nm[k-1] + 1
            Np[k] = Np[k-1]
            j += 1
        k += 1
    while i < n1: # ici, on va rentrer dans UNE SEULE des deux boucles while
        Np[k] = Np[k-1] + 1
        Nm[k] = Nm[k-1]
        k += 1
        i += 1
    while j < n2:
        Nm[k] = Nm[k-1] + 1
        Np[k] = Np[k-1]
        k += 1
        j += 1
    P = Np - Nm
    plt.figure()
    #plt.plot(tau_,Np,'o')
    #plt.plot(tau_,Nm,'+')
    plt.plot(tau_,P)
    #plt.legend(['Np','Nm','P'])
    end_date = time.time()
    print(end_date - start_date)
    plt.show()
    return [tau_,P]

def phi1(x):
    return 0.005*np.exp(-2000*x)

def phi2(x):
    return 1/(x**0.7)


# Export data :

#Hawkes = [simul_Prix(T,1.001,phi1,phi1,1.2) for T in [10,100,500,1000]] 
#d_ = {'tau10':Hawkes[0][0],'data10':Hawkes[0][1]}
#ddf = pd.DataFrame(d_,columns = ['tau10','data10'])
#export_excel = ddf.to_excel(r'C:\\Users\lroussel\Documents\Simulation_Hawkes\export_Hawkes10.xlsx', index = None, header=True)
#d_ = {'tau100':Hawkes[1][0],'data100':Hawkes[1][1]}
#ddf = pd.DataFrame(d_,columns = ['tau100','data100'])
#export_excel = ddf.to_excel(r'C:\\Users\lroussel\Documents\Simulation_Hawkes\export_Hawkes100.xlsx', index = None, header=True)
#d_ = {'tau500':Hawkes[2][0],'data500':Hawkes[2][1]}
#ddf = pd.DataFrame(d_,columns = ['tau500','data500'])
#export_excel = ddf.to_excel(r'C:\\Users\lroussel\Documents\Simulation_Hawkes\export_Hawkes500.xlsx', index = None, header=True)
#d_ = {'tau1000':Hawkes[3][0],'data1000':Hawkes[3][1]}
#ddf = pd.DataFrame(d_,columns = ['tau1000','data1000'])
#export_excel = ddf.to_excel(r'C:\\Users\lroussel\Documents\Simulation_Hawkes\export_Hawkes1000.xlsx', index = None, header=True)
#d_ = {'tau10000':Hawkes[4][0],'data10000':Hawkes[4][1]}
#ddf = pd.DataFrame(d_,columns = ['tau10000','data10000'])
#export_excel = ddf.to_excel(r'C:\\Users\lroussel\Documents\Simulation_Hawkes\export_Hawkes10000.xlsx', index = None, header=True)
#d_ = {'tau15000':Hawkes[5][0],'data15000':Hawkes[5][1]}
#ddf = pd.DataFrame(d_,columns = ['tau15000','data15000'])
#export_excel = ddf.to_excel(r'C:\\Users\lroussel\Documents\Simulation_Hawkes\export_Hawkes15000.xlsx', index = None, header=True)
#d_ = {'tau20000':Hawkes[6][0],'data20000':Hawkes[6][1]}
#ddf = pd.DataFrame(d_,columns = ['tau20000','data20000'])
#export_excel = ddf.to_excel(r'C:\\Users\lroussel\Documents\Simulation_Hawkes\export_Hawkes20000.xlsx', index = None, header=True)



# Ici on va plot la volatilité historique des trajectoires de prix microscopique qu'on a obtenu, afin d'observer une convergence
# vers une vol Rough Heston (trajectoires "rough" avec clusters et retour à la moyenne) + l'effet de levier

def plot_vol_Hawkes(T): 
    os.chdir('C:\\Users\lroussel\Documents\Simulation_Hawkes') 
    data_ = 'export_Hawkes' + str(T)
    df = pd.read_excel(data_ + ".xlsx", index_col=None, header=None)
    df = np.array(df)[1:]
    tau_ = df[:,0]
    x = df[:,1]
    x -= np.min(x)
    x += 1
    y = (x[1:]/x[:-1] - 1)
    #y = x
    sigm = []
    for i in range(9,len(y)): # on prend une fenêtre glissante de 10 sauts
        y_barre = np.mean(y[i-9:i+1])
        sigm.append((1/9)*np.sum((y[i-9:i+1] - y_barre)**2))
    sigm = np.sqrt(sigm)/np.sqrt(T)
    t = tau_[10:]
    plt.figure()
    plt.title("T = " + str(T))
    plt.xlabel("t")
    plt.plot(t,sigm)
    plt.plot(t,x[10:]/np.sqrt(T))
    plt.plot(t,np.mean(sigm)*np.ones(len(t)))
    plt.legend(["Vol","Stock","Vol mean"])
    plt.show()

plot_vol_Hawkes(100)










