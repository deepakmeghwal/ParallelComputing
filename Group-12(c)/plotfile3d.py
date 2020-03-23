import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import csv



mode = int(input('Enter integer [1 -7] for type of DTLZ : '))
obj = int(input('Enter integer [2, 3, 5] for number pf obj functions : '))
trainX = np.genfromtxt('Fun_DTLZ'+str(mode)+'_'+str(obj)+'.txt',delimiter=' ')


X = trainX[:,0]
Y = trainX[:,1]
Z = trainX[:,2]





fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(X, Y, Z, c='r',marker = 'o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.title('DTLZ'+str(mode)+' with 3 objectives 20 variables')
ax.set_xlim3d(0, 1.5)
ax.set_ylim3d(0, 1.5)
ax.set_zlim3d(0, 1.5)




plt.show()
