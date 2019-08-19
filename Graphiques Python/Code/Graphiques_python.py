import numpy as np
import os
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression

os.chdir('C:\\Users\lroussel\Documents\GitHub\Awa_Stages\Data Excel') # Change this


# Extract data :

data_ = 'Facebook_vol'
df = pd.read_excel(data_ + ".xlsx", index_col=None, header=None) # Reading the Excel file
aa = np.flip(np.array(df)[1:])
n = len(aa)


# Study of the Rough feature with realized volatility and looking for the Hurst exponent H :

def m(q,delta): # q-order moment for the log-volatility increments as a function of the lag delta
    diff_log = []
    i = 0
    while i + delta < n:
        diff_log.append(np.log(float(aa[i + delta]/aa[i])))
        i += delta
    res = np.abs(np.array(diff_log))**q
    res = np.mean(res)
    return res

def regr():
    # We plot log(m(q,delta)) as a function of log(delta), for distinct values q
    N = 4 # number of values for q
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
    # We do a linear regression for each q
    x_ = x.reshape(-1,1)
    reg = [LinearRegression().fit(x_,y[i]) for i in range(N)]
    scores = [reg[i].score(x_,y[i]) for i in range(N)]
    print("scores = " + str(scores))
    coefs = [reg[i].coef_ for i in range(N)] # We get the slope for each q
    #print(np.exp(reg[N-1].intercept_/2))   # estimation of the parameter "nu"
    reg_pente = LinearRegression().fit(q_.reshape(-1,1),coefs) # Linear regression for the slope as a function of q
    plt.figure(2)                                              # The slope we will get here corresponds to the Hurst exponent H
    plt.plot(q_,coefs,'o')
    plt.plot(q_,np.dot(q_.reshape(-1,1),reg_pente.coef_) + reg_pente.intercept_)
    plt.xlabel("q")
    plt.ylabel("coef")
    plt.title(data_ + " : pente en fonction de q. On obtient H = " + str(round(float(reg_pente.coef_),2)))
    plt.legend(['coef','linear fit'])
    plt.show()

#regr()



# Get an approximate for H from a FBM sample path :

data_ = 'export_fBMs' 
df = pd.read_excel(data_ + ".xlsx")
aa = np.array(df)[1:]
jj = [aa[:,i] for i in range(9)]
n = len(jj[0])
for i in range(9):
   jj[i] = jj[i] - np.min(jj[i])*2
aa = jj[4] # Change this, putting an integer between 0 and 8
           # 0 for H = 0.1; 1 for H = 0.2; ... ; 8 for H = 0.9
#regr()    



# Study of the SPX ATM skew :

def regr_skew():
    os.chdir('C:\\Users\lroussel\Documents\GitHub\Awa_Stages\Data Excel') # Change this
    data_ = 'LiftedHestonVolSurface' 
    df = pd.read_excel(data_ + ".xlsx")
    df = np.array(df)
    values = df[1:,1:]
    maturities = df[1:,0]
    logMoneyness = df[0,1:]
    n1 = len(maturities)
    n2 = len(logMoneyness)
    skew = np.zeros((n1,n2-1))
    for i in range(n2-1):
        skew[:,i] = (values[:,i+1] - values[:,i])/(logMoneyness[i+1] - logMoneyness[i])
    skew = skew[:,9] # get ATM skew
    logMaturity = np.log(maturities)[:8]
    logSkew = np.log(np.abs(skew))[:8]
    reg_skew = LinearRegression().fit(logMaturity.reshape(-1,1),logSkew)
    plt.figure(1)
    plt.plot(np.exp(logMaturity),np.exp(logSkew),'o')
    plt.plot(np.exp(logMaturity),np.exp(np.dot(logMaturity.reshape(-1,1),reg_skew.coef_) + reg_skew.intercept_))
    plt.title("SPX : ATM skew en fonction de la maturité. alpha = " + str(round(-float(reg_skew.coef_),2)))
    plt.legend(['ATM skew','fit'])
    plt.show()

#regr_skew()

