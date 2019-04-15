import numpy as np
from scipy.special import gamma
from matplotlib import pyplot as plt

def BM(borne_min,T,n):
    dt = (T-borne_min)/n
    N = np.random.normal(0,1,n)
    W = np.zeros(n+1-int(T/dt))
    W2 = np.zeros(int(T/dt))
    for i in range(1,n+1-int(T/dt)):
        W[i] = W[i-1] + np.sqrt(dt)*N[i-1]
    W = np.flip(W)
    W2[0] = W[n-int(T/dt)] + np.sqrt(dt)*N[n-int(T/dt)]
    for i in range(n+1-int(T/dt),n):
        W2[i-n+int(T/dt)] = W2[i-n+int(T/dt)-1] + np.sqrt(dt)*N[i]
    W = np.concatenate((W,W2))
    return W

def f1(s,t,gam):
    return 1/((t - s)**gam)

def f2(s,gam):
    return f1(s,0,gam)

def fBM(H,borne_min,T,n): # représentation de Mandelbrot - Van Ness
    dt = (T-borne_min)/n
    gam = 1/2 - H
    C = np.sqrt(2*H*gamma(3/2 - H)/(gamma(H + 1/2)*gamma(2 - 2*H)))
    W = BM(borne_min,T,n)
    fW = np.zeros(int(T/dt)+1) 
    C2 = f2(borne_min,gam)*(W[1] - W[0])
    for i in range(1,n-int(T/dt)):
        C2 += f2(borne_min + i*dt,gam)*(W[i+1] - W[i])
    fW[0] = C2
    for i in range(1,int(T/dt)+1):
        fW[i] = fW[i-1] + f1(borne_min + (n-int(T/dt)+i-1)*dt,borne_min + (n-int(T/dt)+i)*dt,gam)*(W[n-int(T/dt)+i] - W[n-int(T/dt)+i-1])
    fW = (fW - C2)
    return fW

def plot_fBM():
    T = 2
    borne_min = -198
    n = 100000
    dt = (T-borne_min)/n
    H = [0.1, 0.2, 0.5, 0.7]
    fBMs = [fBM(h,borne_min,T,n) for h in H]
    X = [borne_min + (n-int(T/dt)+i)*dt for i in range(int(T/dt)+1)]
    plt.figure()
    for m in fBMs:
        plt.plot(X,m)
    plt.xlabel("t")
    plt.ylabel("W(t,H)")
    plt.title("Simulation fBM par représentation de Mandelbrot - Van Ness")
    plt.legend(["H = " + str(h) for h in H])
    plt.show()

plot_fBM()


