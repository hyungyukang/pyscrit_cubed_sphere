import numpy as np
#-----------------------------------------------------------------------
#
#  Sphere to Cartesian
#
#-----------------------------------------------------------------------
def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

#-----------------------------------------------------------------------
#
#  Legendre-Gauss-Lobatto (LGL) points and weights
#
#  input : N = # of gridpoints
#  output: xg,wg
#          xg : gaussian points
#          wg : gaussian weights at each gaussian point
#
#-----------------------------------------------------------------------

def legendre_gauss_lobatto(a,b,N):
    apb= (a+b)/2.0
    bma= (b-a)/2.0
    
    #nh= np.floor(N+1)/2.0
    #Np= N+1
    
    xg= np.zeros([N])
    wg= np.zeros([N])
    xx=0
    
    xg[0]= -1
    wg[0]=2.0/(N*(N-1))
    xg[N-1]= +1.0
    wg[N-1]=2.0/(N*(N-1))

    for i in range(1,N-1):
        #x= cos( (2*i-1)*pi/(2*N+1) );
        x=(1.0-(3.0*(N-2.0))/(8.0*(N-1)**3)) * np.cos(np.pi*(4.0*(i+1)-3)/(4.0*(N-1)+1))
    
        for j in range(100):
            P0=1
            P1=0
            dp=1
            dpp=1
            dppp=1
            for m in range(N):
                n = m+1  
                P2=P1
                P1=P0
                dp1= dp
                dpp1= dpp
                dppp1= dppp
         
                P0= (2.0*n-1)*x*P1/n - (n-1)*P2/n
                dp= n*(P1-x*P0)/(1.0-x*x) # the same with pp= (n/(x*x-1))*(x*P0-P1)
                dpp= ( 2.0*x*dp - n*(n+1)*P0   )/(1.0-x*x)
                dppp= ( 2.0*x*dpp-(n*(n+1)-2.0)*dp)/(1.0-x*x)
            xx= x-(2.0*dp1*dpp1)/(2.0*dpp1*dpp1 - dp1*dppp1)
            if (np.abs(xx-x)<10.e-15):
               break
            x= xx
 
        xg[N-i]= x*bma+apb
        xg[  i]= -xg[N-i]
        wg[N-i]= 2.0/(N*(N-1)*P1*P1)
        wg[  i]=  wg[N-i]

    #x=xg
    #w=wg

    return xg,wg

#-----------------------------------------------------------------------
#
#  Cubed sphere grid
#
#-----------------------------------------------------------------------
def function_Cubed_sphere_grid_equiang (No,Ne,IB):
   Ng = No + 1  # GLL points
   Nge = Ng*Ne
   Nfaces = 6
   
   de_x = 2 / (Ne)  # delta element (-1 ~ 1)
   de_ang = (np.pi/2) / (Ne) # delta element (-pi/4 ~ pi/4)
   
   dx = 2/(IB-1)        # delta x (-1 ~ +1)
   dthe = (np.pi/2)/(IB-1)  # delta pi (-pi/4 ~ pi/4)
   a = 1.0  # Unit Earth
   dist_limit = 1.0e-10 
   pi2 = np.pi*2.0
   

   xg= np.zeros(Ng)
   wg= np.zeros(Ng)
   
   #%%%%%%%%%%%%   Call GLL subroutine %%%%%%%%%%%%%
   xg, wg = legendre_gauss_lobatto(-1,+1,Ng)
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ELEp = np.zeros(Ne)
   Eratio = np.zeros([Ng,Ne])
   xgrid = np.zeros(IB)
   xcoord = np.zeros([Nge,IB])
   cart_y = np.zeros([Nge])
   cart_x = np.zeros([IB,Nge])
   cart_r = np.zeros([IB,Nge])
   sphere_lon = np.zeros([IB,Nge,Nfaces])
   sphere_lat = np.zeros([IB,Nge,Nfaces])
   x = np.zeros([IB,Nge,Nfaces])
   y = np.zeros([IB,Nge,Nfaces])
   z = np.zeros([IB,Nge,Nfaces])
   
   # For elements line  ( -pi/4 ~ pi/4 )
   for n in range (Ne):
         ELEp[n] = -np.pi/4 + n*de_ang
   #     ELEp[n] = -1 + n*de_x # For equidis
   
   for k in range (Ng):
      for n in range (Ne):
          Eratio[k,n] = ELEp[n]+de_ang/2  + xg[k]*(de_ang/2)
      #   Eratio[n,k] = ELEp[n]+de_x/2  + xg[k]*(de_x/2) # For equidis
   
   for i in range(IB):
        xgrid[i] = -np.pi/4 + i*dthe 
   #    xgrid[i] = -1.0 + i*dx  # For equidis
   
   ic = 0
   for n in range(Ne):
       for k in range(Ng):
           for i in range(IB):
               xcoord[ic,i] = Eratio[k,n]
           ic = ic + 1

   ic = 0
   for n in range(Ne):
       for k in range(Ng):
           cart_y[ic] = a*np.tan(Eratio[k,n])
           ic = ic + 1
   
   for j in range(Nge):
       for i in range(IB):
           cart_x[i,j]=a*np.tan(xgrid[i]) 
           cart_r[i,j] = np.sqrt(a**2.0+ cart_x[i,j]**2.0+cart_y[j]**2.0);
   
   for ip in range(Nfaces):
       for j in range(Nge):
           for i in range(IB):
               if (ip == 0):
                  sphere_lat[i,j,ip]= np.arcsin(cart_y[j]/cart_r[i,j])
                  sphere_lon[i,j,ip]= np.arctan2(cart_x[i,j], a)
               elif(ip == 1):
                  sphere_lat[i,j,ip]= np.arcsin(cart_y[j]/cart_r[i,j])
                  sphere_lon[i,j,ip]= np.arctan2(a,-cart_x[i,j])
               elif (ip == 2):
                  sphere_lat[i,j,ip]= np.arcsin(cart_y[j]/cart_r[i,j])
                  sphere_lon[i,j,ip]= np.arctan2(-cart_x[i,j],-a)
               elif (ip == 3):
                  sphere_lat[i,j,ip]= np.arcsin(cart_y[j]/cart_r[i,j])
                  sphere_lon[i,j,ip]= np.arctan2(-a,cart_x[i,j])
               elif (ip == 4):
                  yy = np.abs(cart_y[j])
                  xx = np.abs(cart_x[i,j])
                  if (xx > dist_limit or yy > dist_limit):
                     sphere_lon[i,j,ip]= np.arctan2(cart_x[i,j],-cart_y[j])
                  else:
                     sphere_lon[i,j,ip]= 0.0
                  sphere_lat[i,j,ip]=  np.arcsin( a/cart_r[i,j]);
               elif (ip == 5):
                  yy = np.abs(cart_y[j]);
                  xx = np.abs(cart_x[i,j]);
                  if (xx > dist_limit or yy > dist_limit): 
                     sphere_lon[i,j,ip]= np.arctan2(cart_x[i,j], cart_y[j])
                  else:
                     sphere_lon[i,j,ip]= 0.0
                  sphere_lat[i,j,ip]=  np.arcsin(-a/cart_r[i,j])
   
               if (sphere_lon[i,j,ip] < 0.0 ):
                   sphere_lon[i,j,ip] = sphere_lon[i,j,ip] + pi2
   
   
   for ip in range(Nfaces):
       for k in range(Nge):
           for i in range(IB):
               x[i,k,ip],y[i,k,ip],z[i,k,ip] = sph2cart(sphere_lon[i,k,ip],sphere_lat[i,k,ip],1.0) 
 
   return x,y,z
