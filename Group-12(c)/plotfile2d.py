import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import csv



mode = int(input('Enter integer [1 -7] for type of DTLZ : '))
obj = int(input('Enter integer [2, 3, 5] for number pf obj functions : '))
trainX = np.genfromtxt('Fun_DTLZ'+str(mode)+'_'+str(obj)+'.txt',delimiter=' ')


X = trainX[:,0]
Y = trainX[:,1]



fig = plt.figure()
plt.scatter(X, Y, c= 'red', alpha=0.1)
plt.xlabel('Function 1')
plt.ylabel('Function 2')
plt.xlim(0,70)
plt.ylim(0,70)
plt.title('DTLZ'+str(mode)+' with 2 objectives 20 variables')
plt.show()
