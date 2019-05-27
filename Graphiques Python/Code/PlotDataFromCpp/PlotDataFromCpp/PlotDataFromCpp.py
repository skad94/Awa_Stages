import numpy as np
import os
from matplotlib import pyplot as plt

os.chdir('C:\\Users\\lroussel\\Documents\\GitHub\\Awa_Stages\\Library\\repos\\Math') # Change this

# Plotting a FBM sample path which was generated using C++ :

dataCholesky03 = 'SimulationFBM_Cholesky_03.txt' 
dataCholesky08 = 'SimulationFBM_Cholesky_08.txt'
dataKL03 = 'SimulationFBM_KL_03.txt'
dataKL08 = 'SimulationFBM_KL_08.txt'
f = open(dataKL03, 'r+') # replace the first argument with one of the 4 names above
FBM_path = f.read().split(',')[:-1]
f.close()
FBM_path = np.array(FBM_path).astype(np.float)
plt.figure()
plt.plot(FBM_path)
plt.show()
