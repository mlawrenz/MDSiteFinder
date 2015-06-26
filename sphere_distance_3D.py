from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy
import scipy.spatial as sp

x=numpy.linspace(0,20,20)
y=numpy.linspace(0,20,20)
z=numpy.linspace(0,20,20)
a,b,c=numpy.meshgrid(x,y,z)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(a,b,c, alpha=0.2, c='k')
cutoff=2
x=x[abs(x-4)<cutoff]
y=y[abs(y-4)<cutoff]
z=z[abs(z-4)<cutoff]
X,Y,Z=numpy.meshgrid(x,y,z)
ax.scatter(X,Y,Z, c='r', alpha=0.5)
data=numpy.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
print data.shape
center=numpy.array([4,4,4])
distance=sp.distance.cdist(data, center.reshape(1,-1)).ravel()
points_in_sphere=data[distance < cutoff]
print points_in_sphere

import pdb
pdb.set_trace()
for (i,j,k) in points_in_sphere:
    ax.scatter(i,j,k, c='r')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

