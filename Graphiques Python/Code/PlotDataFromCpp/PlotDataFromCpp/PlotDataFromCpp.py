import numpy as np
import os
from matplotlib import pyplot as plt


os.chdir('C:\\Users\\lroussel\\Documents\\GitHub\\Awa_Stages\\Library\\repos\\Math') # Change this

# Plotting a FBM sample path which was generated using C++ :

dataLiftedHeston = 'LiftedHestonPrice.txt' 

fLiftedHeston = open(dataLiftedHeston, 'r+')
pathLiftedHeston = fLiftedHeston.read().split(',')[:-1]
fLiftedHeston.close()
pathLiftedHeston = np.array(pathLiftedHeston).astype(np.float)
plt.figure()
plt.plot(pathLiftedHeston)
plt.title("Trajectoire du sous-jacent pour le modèle Lifted Heston")
plt.xlabel('t')
plt.ylabel('S')
#plt.legend(['rHeston','rBergomi'])
plt.show()

