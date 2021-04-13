#from mpl_toolkits import mplot3d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from function_Cubed_sphere_grid_equiang import function_Cubed_sphere_grid_equiang

#------------------------------------------------------------------------
#
#  Python script for the cubed sphere grid (translated from MATLAB codes)
#
#                      H.-G. Kang (APR 13 2021)
#
#------------------------------------------------------------------------


No =3    # Polynomial Order
Ne = 4       # #Element
IB = 501  # Resolution

Ng = No + 1  # GLL points
Nge = Ng*Ne
Nfaces = 6

#------------------------------------------------------------------------

x = np.zeros([IB,Nge,Nfaces])
y = np.zeros([IB,Nge,Nfaces])
z = np.zeros([IB,Nge,Nfaces])

x,y,z = function_Cubed_sphere_grid_equiang (No,Ne,IB)

fig,ax = plt.subplots(1,figsize=(6,5),constrained_layout=True)
#fig = plt.subplots(figsize=(6,5),constrained_layout=True)
#------------------------------------------------------------------------

# Display ----------------------------------------------------

#ax = plt.axes(projection="3d")

for ip in range(Nfaces):
    for k in range(Nge):
        #ax.scatter3D(x,y,z,s=10)
        ax.scatter(x[:,k,ip],y[:,k,ip],z[:,k,ip],c='k')
        ax.scatter(y[:,k,ip],x[:,k,ip],z[:,k,ip],c='k')
        ax.scatter(z[:,k,ip],x[:,k,ip],y[:,k,ip],c='k')

#------------------------------------------------------------------------
# Another set for red colored grid

No = 1    # Polynomial Order
#Ne = 3       # #Element
IB = 1001  # Resolution

Ng = No + 1  # GLL points
Nge = Ng*Ne

xr = np.zeros([IB,Nge,Nfaces])
yr = np.zeros([IB,Nge,Nfaces])
zr = np.zeros([IB,Nge,Nfaces])

xr,yr,zr = function_Cubed_sphere_grid_equiang (No,Ne,IB)


for ip in range(Nfaces):
    for k in range(Nge):
        ax.scatter(xr[:,k,ip],yr[:,k,ip],zr[:,k,ip],color='r')
        ax.scatter(yr[:,k,ip],xr[:,k,ip],zr[:,k,ip],color='r')
        ax.scatter(zr[:,k,ip],xr[:,k,ip],yr[:,k,ip],color='r')

ax.axis('equal')
plt.savefig("fig_test.png")
