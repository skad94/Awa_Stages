import numpy as np
from scipy.special import gamma
from matplotlib import pyplot as plt
import scipy.integrate as integrate
import time

start_time = time.time()


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
    elapsed_time = time.time() - start_time
    print("time = " + str(elapsed_time))
    plt.show()


def my_cos(t,k,T):
    return np.cos(np.pi*t*k/T)

def my_sin(t,k,T):
    return np.sin(np.pi*t*k/T)

def c0_KL(H,T): # on définit les coefficients qui interviennent dans la méthode par "series expansion"
    if H < 1/2:
        return 0
    return H*(T**(2*H-2))

def ck_KL(H,T,k):
    def f_(u):
        return (u**(2*H))*my_cos(u,k,T)
    def g_(u):
        return (u**(2*H-2))*my_cos(u,k,T)
    if H < 1/2:
        return (2/T)*integrate.quad(f_, 0, T)[0]
    return (-4*H*(2*H-1)*T/((k*np.pi)**2))*integrate.quad(g_, 0, T)[0]

def fBM3(H,T,Nmax,n): # simulation des trajectoires par "series expansion"
    dt = T/n
    fW = np.zeros(n+1)
    N = np.random.normal(0,1,2*Nmax+1)
    ck = [ck_KL(H,T,k) for k in range(1,Nmax+1)]
    for i in range(1,n+1):
        res = np.sqrt(c0_KL(H,T))*i*dt*N[0]
        for k in range(1,Nmax+1):
            tmp = np.sqrt(-ck[k-1]/2)
            res += tmp*(my_sin(i*dt,k,T)*N[2*k] + (1 - my_cos(i*dt,k,T))*N[2*k-1])
        fW[i] = res
    return fW
        
def plot_fBM3(H,T,Nmax,n):
    dt = T/n
    fBMs = [fBM3(h,T,Nmax,n) for h in H]
    X = [i*dt for i in range(n+1)]
    plt.figure()
    for m in fBMs:
        plt.plot(X,m)
    plt.xlabel("t")
    plt.ylabel("W(t,H)")
    plt.title("Simulation fBM par \"series expansion\"")
    plt.legend(["H = " + str(h) for h in H])
    plt.show()



#plot_fBM([0.3,0.5,0.7],-198,2,100000)
#plot_fBM2([0.3,0.5,0.7],2,1000)
#plot_fBM3([0.3,0.55,0.7],1,240,300)




















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


#def trapezes_stochastique(a,b,K,f,W): # calcul d'une intégrale stochastique par la méthode des trapèzes
#    delta = (b - a)/K
#    res = 0
#    for i in range(K):
#        res += (f(a + (i+1)*delta) + f(a + i*delta))*(W[a + (i+1)*delta] - W[a + i*delta])/2
#    return res

#def trapezes(T,dt,f):
#    res = 0.9*dt*(f(0.1*dt) + f(dt))/2
#    for i in range(1,int(T/dt)):
#        res += dt*(f((i+1)*dt) + f(i*dt))/2
#    return res

