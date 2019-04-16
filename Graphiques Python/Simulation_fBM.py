import numpy as np
from scipy.special import gamma
from matplotlib import pyplot as plt

def BM(borne_min,T,n):
    dt = (T-borne_min)/n
    N = np.random.normal(0,1,n)
    W = np.zeros(n+1)
    for i in range(1,n+1):
        W[i] = W[i-1] + np.sqrt(dt)*N[i-1]
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
    W -= W[n-int(T/dt)] # on fait en sorte que le processus soit égal à 0 en 0
    fW = np.zeros(int(T/dt)+1) 
    C2 = f2(borne_min,gam)*(W[1] - W[0])
    for i in range(1,n-int(T/dt)-1):
        C2 += (1/2)*(f2(borne_min + i*dt,gam)+f2(borne_min + (i+1)*dt,gam))*(W[i+1] - W[i])
    C2 += f2(borne_min + (n-int(T/dt)-1)*dt,gam)*(W[n-int(T/dt)] - W[n-int(T/dt)-1])
    fW[0] = C2
    for i in range(1,int(T/dt)+1):
        fW[i] = fW[i-1] + f1(dt/2,dt,gam)*(W[n-int(T/dt)+i] - W[n-int(T/dt)+i-1])
    fW = C*(fW - C2)
    return fW

def plot_fBM(H,borne_min,T,n):
    dt = (T-borne_min)/n
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

def cov_fBM(s,t,H): # fonction de covariance du fBM
    return (s**(2*H) + t**(2*H) - np.abs(t-s)**(2*H))/2

def fBM2(H,T,n): # décomposition de Cholesky
    dt = T/n
    N = np.random.normal(0,1,n)
    COV_ = np.zeros((n,n))
    for i in range(1,n+1):
        for j in range(1,n+1):
            COV_[i-1,j-1] = cov_fBM(i*dt,j*dt,H)
    SIGM_ = np.linalg.cholesky(COV_)
    N = np.concatenate(([0],np.dot(SIGM_,N)))
    return N

def plot_fBM2(H,T,n):
    dt = T/n
    fBMs = [fBM2(h,T,n) for h in H]
    X = [i*dt for i in range(n+1)]
    plt.figure()
    for m in fBMs:
        plt.plot(X,m)
    plt.xlabel("t")
    plt.ylabel("W(t,H)")
    plt.title("Simulation fBM par décomposition de Cholesky")
    plt.legend(["H = " + str(h) for h in H])
    plt.show()

#plot_fBM([0.1,0.3,0.5,0.8],-198,2,100000)
#plot_fBM2([0.2,0.3,0.5,0.8],2,1000)






#def verif(): # Fonction pour vérifier que le processus a la bonne covariance (version non définitive...)
#    H = [0.7]
#    M = 10
#    borne_min = -98
#    T = 2
#    n = 500
#    dt1 = (T-borne_min)/n
#    dt2 = T/n
#    X1 = [borne_min + (n-int(T/dt1)+i)*dt1 for i in range(int(T/dt1)+1)]
#    X2 = [i*dt2 for i in range(n+1)]
#    c = cov_fBM(X2[2],X2[90],H[0])
#    res = 0
#    for i in range(M):
#        tmp = fBM2(H[0],T,n)
#        res += tmp[2]*tmp[90]
#    res /= M
#    print(np.abs(res - c))

#verif()


