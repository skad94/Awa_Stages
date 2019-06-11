import numpy as np
import os
from matplotlib import pyplot as plt

os.chdir('C:\\Users\\lroussel\\Documents\\GitHub\\Awa_Stages\\Library\\repos\\Math') # Change this

# Plotting a FBM sample path which was generated using C++ :

dataCholesky03 = 'SimulationFBM_Cholesky_03.txt' 
dataCholesky08 = 'SimulationFBM_Cholesky_08.txt'
dataKL03 = 'SimulationFBM_KL_03.txt'
dataKL08 = 'SimulationFBM_KL_08.txt'
data_rHeston = 'rHestonPrice.txt'
data_rBergomi = 'rBergomiPrice.txt'
fHeston = open(data_rHeston, 'r+')
fBergomi = open(data_rBergomi, 'r+')
path_Heston = fHeston.read().split(',')[:-1]
path_Bergomi = fBergomi.read().split(',')[:-1]
fHeston.close()
fBergomi.close()
path_Heston = np.array(path_Heston).astype(np.float)
path_Bergomi = np.array(path_Bergomi).astype(np.float)
plt.figure()
plt.plot(path_Heston)
plt.plot(path_Bergomi)
plt.title("Trajectoire du sous-jacent pour les modèles rHeston et rBergomi")
plt.xlabel('t')
plt.ylabel('S')
plt.legend(['rHeston','rBergomi'])
plt.show()
