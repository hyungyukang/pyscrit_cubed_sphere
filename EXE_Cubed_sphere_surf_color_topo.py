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
Ng = No + 1  # GLL points
Ne = 4       # #Element
Nge = Ng*Ne
IB = 501  # Resolution
Nfaces = 6

pi180  = np.pi/180

x = np.zeros([IB,Nge,Nfaces])
y = np.zeros([IB,Nge,Nfaces])
z = np.zeros([IB,Nge,Nfaces])

x,y,z = function_Cubed_sphere_grid_equiang (No,Ne,IB)


# Display ----------------------------------------------------

fig,ax = plt.subplots(1,figsize=(6,5),constrained_layout=True)

for ip in range(Nfaces):
    for k in range(Nge):
        ax.scatter(x[:,k,ip],y[:,k,ip],z[:,k,ip],color='k')
        ax.scatter(y[:,k,ip],x[:,k,ip],z[:,k,ip],color='k')
        ax.scatter(z[:,k,ip],x[:,k,ip],y[:,k,ip],color='k')

ax.axis('equal')
plt.savefig("fig_test.png")
