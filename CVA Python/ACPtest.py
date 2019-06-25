import pandas as pd
import numpy as np
import xlrd
from xlwt import Workbook, Formula
import matplotlib.pyplot as plt
from math import *
from scipy import interpolate
from sklearn.linear_model import LinearRegression
import scipy.integrate as integrate
import scipy.special as special
from random import gauss
import time
import csv

def invert(x):
    y = x 
    n = len(x)
    for i in range(n):
        x[i]=y[n-i-1]
    return x

def InterpolationLinear(x , map, convention = 'linear'):
   if (x==0):
        return 0
   time = map.keys()
   n = len(time)
   r2  = 0
   if (x < float(list(time)[0])):
       return map[list(time)[0]]
   if (x >= float(list(time)[n-1])):
       return map[list(time)[n-1]]
   for i in range(n):
       t2 = list(time)[i]
       r2 = map[t2]
       if ( x < float(t2)):
           t1 = list(time)[i-1]
           r1 = map[t1]
           return (x- t2) / (t1 - t2) * r1 + (x - t1) / (t2 - t1) * r2
       
def Interpolation2(x, PD):
    n = PD.shape[1]
    r2 = 0
    if (x <= PD[0][0]):
        return PD[1][0]
    if (x>= PD[0][PD.shape[1]-1]):
        return PD[1][PD.shape[1]-1]
    for i in range(n):
        t2=PD[0][i]
        r2 = PD[1][i]
        if (t2 >= x ):
            t1 = PD[0][i-1]
            r1 = PD[1][i-1]
            return (x - t2) / (t1 - t2) * r1 + (x - t1) / (t2 - t1) * r2


def Rsquared(x,y,f):
    mean = 0
    for i in range(x.shape[0]):
        mean += y[i]
    mean =mean / x.shape[0]
    denom = 0
    num = 0
    for i in range(x.shape[0]):
        num += pow(f(x[i]) - y[i],2 )
        denom += pow(y[i] - mean,2)
    return (1 - num/denom)

        
def IntegrateRectangle(f,a,b,N):
    h = (b-a)/N
    inte=0
    for i in range(0,N):
        x_i = a+ i*h
        x_imore1=a+(i+1)*h
        inte += f( (x_i+x_imore1)*0.5)
    return h*inte

    
def Calibrage(name):#Return tupple (time,fit[1,2,3],information rate,RSquared[1,2,3])
    document = xlrd.open_workbook("C:/Users/alepeltier/Documents/CVA/"+name)
    feuille_1 = document.sheet_by_index(0)
    cols = feuille_1.ncols
    rows = feuille_1.nrows
    time=np.zeros(cols-1)
    for i in range(0,time.shape[0]):
        time[i]=feuille_1.cell_value(0,i+1)
    time[0]=0
    Data = np.eye(rows-1,cols-1)
    for i in range(1, rows):
        for j in range(1, cols):
            Data[i-1,j-1] = feuille_1.cell_value(i,j)
    f_0 = Data[Data.shape[0]-1,:]/100
    Diff=np.eye(rows-2,cols-1)
    for i in range(0,rows-2):
        for j in range(0,cols-1):
            Diff[i,j]=Data[i+1,j]-Data[i,j]
            
    
    Expectation  = np.zeros(Diff.shape[1])   
    for j in range(0,Diff.shape[1]):#Colonne
        for i in range(0,Diff.shape[0]):#Ligne
             Expectation[j] += Diff[i,j]
        Expectation[j] = Expectation[j] / (Diff.shape[0]-1)
    
    CovMat = np.eye(Diff.shape[1],Diff.shape[1])
    for i in range(0,Diff.shape[1]):
        for j in range(0,Diff.shape[1]):
            CovMat[i,j]=0
            for k in range(0,Diff.shape[0]):
                CovMat[i][j] += (Diff[k][i] - Expectation[i])*(Diff[k][j] - Expectation[j])
            CovMat[i][j]=CovMat[i][j] / (Diff.shape[0]) * (252/10000)    
    A=np.linalg.eig(CovMat)    
    cumulative = np.zeros(A[0].shape)
    sum=0.0
    for i in range(0,cumulative.shape[0]):
        for j in range(0,i+1):
            cumulative[i]+= (A[0])[j]
        cumulative[i]=cumulative[i]/np.sum(A[0])  
    v1 = sqrt(A[0][0])*(A[1])[:,0]
    v2 = sqrt(A[0][1])*(A[1])[:,1]
    v3 = -sqrt(A[0][2])*(A[1])[:,2]
        
    fit1 = np.polyfit(time,v1,3)
    fit2 = np.polyfit(time,v2,3)
    fit3 = np.polyfit(time,v3,3)
#    fit1 = [0.0064,0,0,0]
#    fit2 = [-0.0036, -0.00057 , 0.00012 , -0.0000036]
#    fit3 = [-0.0048 , 0.0018 , -0.00014 , 0.0000031]
#    fit1 = [0,0,0,0.0064]
#    fit2 = [--0.0000036, 0.00012 , -0.00057 , -0.0036 ]
#    fit3 = [ 0.0000031 , -0.00014 , 0.0018 , -0.0048]
    fit=(fit1,fit2,fit3)
#    erreur1 = Rsquared(time, v1,f1)    
#    erreur2 = Rsquared(time, v2,f2)
#    erreur3 = Rsquared(time, v3,f3)
#    erreur=(erreur1,erreur2,erreur3)
    
    def f1(x):
        return fit[0][3]+fit[0][2]*x+fit[0][1]*pow(x,2)+fit[0][0]*pow(x,3)

    def f2(x):
        return fit[1][3]+fit[1][2]*x+fit[1][1]*pow(x,2)+fit[1][0]*pow(x,3)

    def f3(x):
        return fit[2][3]+fit[2][2]*x+fit[2][1]*pow(x,2)+fit[2][0]*pow(x,3)
    
    def m_barre(Tau,IntegrateMethod=IntegrateRectangle):
        N=10
        if Tau==0:
            return 0
        else:
            M1 = f1(Tau)* IntegrateMethod(f1,0,Tau,N)
            M2 = f2(Tau)* IntegrateMethod(f2,0,Tau,N)
            M3 = f3(Tau)* IntegrateMethod(f3,0,Tau,N)
            return M1+M2+M3
    
    m=np.zeros(time.shape[0])    
    for i in range(time.shape[0]):
        m[i]=m_barre(time[i])
        
    f = np.zeros((3,time.shape[0]))
    for i in range(time.shape[0]):
        f[0][i]=f1(time[i])
        f[1][i]=f2(time[i])
        f[2][i]=f3(time[i])
            
    calib=(time,f_0,fit,cumulative[3],m,f)     
    return calib
    

def Euler(calib,T,dTau = 0.01):
    def f1(x):
        return calib[2][0][3]+calib[2][0][2]*x+calib[2][0][1]*pow(x,2)+calib[2][0][0]*pow(x,3)

    def f2(x):
        return calib[2][1][3]+calib[2][1][2]*x+calib[2][1][1]*pow(x,2)+calib[2][1][0]*pow(x,3)

    def f3(x):
        return calib[2][2][3]+calib[2][2][2]*x+calib[2][2][1]*pow(x,2)+calib[2][2][0]*pow(x,3)
    
#    
#    def m_barre(Tau,dTau,IntegrateMethod=IntegrateRectangle):
#        N=10
#        if Tau==0:
#            return 0
#        else:
#            M1 = f1(Tau)* IntegrateMethod(f1,0,Tau,N)
#            M2 = f2(Tau)* IntegrateMethod(f2,0,Tau,N)
#            M3 = f3(Tau)* IntegrateMethod(f3,0,Tau,N)
#            return M1+M2+M3
##        
    time = calib[0]
    N = int(T/dTau)
    Traj = np.zeros((N,time.shape[0]))
    Traj[0,:]=calib[1]
    for i in range(1,N):
        n1 = gauss(0,1)
        n2 = gauss(0,1)
        n3 = gauss(0,1)
#        n1 = np.random.randn(1)
#        n2 = np.random.randn(1)
#        n3 = np.random.randn(1)
        for j in range(0,time.shape[0]):
            Traj[i,j] = Traj[i-1,j]
            Traj[i,j] += sqrt(dTau) *(n1*calib[5][0][j] + n2*calib[5][1][j] + n3*calib[5][2][j])
            if j == (time.shape[0]-1):
                Traj[i,j] += dTau * (calib[4][j] + ((Traj[i-1,j]-Traj[i-1,j-1]) / (time[j]-time[j-1])))
            else:
                Traj[i,j] += dTau * (calib[4][j] + ((Traj[i-1,j+1]-Traj[i-1,j]) / (time[j+1]-time[j])))
    return Traj  

def ZC_HJM(t,T,calib, M=100,dTau = 0.01):
#    zc_conv = np.zeros(M-1)#tableau for the convergence
    if (t==T):
        return 1
    zero_coupon = 0
    for j in range(M):
#        if(j%10==0 and j>0):
#            print(j,"   ",zero_coupon/j)
        Traj = Euler(calib,T,dTau)
        sum = 0
        for i in range(int(t/dTau),int(T/dTau)):
            sum+=Traj[i,0]
        zero_coupon += exp (- dTau* sum)
#        if (j>0):
#            zc_conv[j-1] = (zero_coupon/(j+1))
    return zero_coupon


def Swap(t,T,frequency, nominal, fix_rate,calib,M):
    k=0
    while (k < t):
        k+=frequency
    if (k==0):
        k += frequency
    float_leg= ZC_HJM(t,k,calib,M) - ZC_HJM(t,T,calib,M)
    fix_leg=0
    for i in np.arange(k,T+frequency,frequency):
        fix_leg+= ZC_HJM(t,i,calib,M) * frequency * fix_rate
    return nominal * (fix_leg - float_leg)

def curve (fichier):
    with open('C:/users/alepeltier/Documents/GitHub/Awa_Stages/Data/'+fichier+'.csv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        i = 0 
        for row in spamreader:
            if (i == 0):
                map = { row[0] : row[1] }
            else:
                map[float(row[0])] = float(row[1]);
            i=i+1
    return map

def Actualisation (curve,t):
    if (t==0):
        return 1
    else:
        return pow(1/(1+InterpolationLinear(t,curve)),t)
    
def CVA(start_date , maturity , frequency , nominal, fixe_rate ,N,M,curve):
    nom ="ACP_Antoine.xlsx"
    calib=Calibrage(nom)
    cva = 0
    delta_t = maturity/N
    for i in range (N+1):
        print("durée", delta_t*i)
        prix = Swap( delta_t * i, maturity , frequency , nominal , fixe_rate , calib , M )
        print("prix", prix)
        print("proba defaut",Interpolation2(delta_t *(i+1),PD) - Interpolation2(delta_t * i,PD))
        print("acutalisation",Actualisation(curve,delta_t*i))
        #print(Actualisation(curve,delta_t*i))
        cva += Actualisation(curve,delta_t*i) * prix * (Interpolation2(delta_t *(i+1),PD) - Interpolation2(delta_t * i,PD))
    return cva
    
    

nom ="ACP_Antoine.xlsx"
nom2="HJM Model - PCA.xlsm"
calib=Calibrage(nom)

#zc=ZC_HJM(0,2,calib,0.01,1000)
#print(zc)


#Traj = Euler(calib,10,0.01)

#prix = Swap(0 , 8 , 1 , 100 , - 0.2 , calib)
#print(prix)
  

discount = curve('Output') 

PD = np.array([[0,0.5,1,2,3,4,5,7,10,20,30],[1,0.99,0.99,0.98,0.96,0.94,0.93,0.90,0.86,0.74,0.63]])
PD[1] = 1 - PD[1]

#print(Interpolation2(7.2,PD))


fixe_rate = - 0.2/100
start_date= 0
maturity = 5
frequency = 1
nominal = 100



cva = CVA(start_date , maturity , frequency , nominal, fixe_rate , 2 , 2, discount)
print("cva",cva)
prix = Swap( start_date, maturity , frequency , nominal , fixe_rate , calib , 50 )
print("prix",prix)




    
