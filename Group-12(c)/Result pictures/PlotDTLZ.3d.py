import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import csv

trainX = np.genfromtxt('DTLZ7.3D.pf',delimiter=' ')
X = trainX[:,0]
Y = trainX[:,1]
Z = trainX[:,2]




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(X, Y, Z, c='r',marker = 'o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
