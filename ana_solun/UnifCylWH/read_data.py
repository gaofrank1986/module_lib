import numpy as np
import bessl_func as bf
from sympy import mpmath as mp
import os,sys
f = open("./DPOINT.txt")
(h,wk) = f.readline().split()
h = float(h)
wk = float(wk)

num_pts = int(f.readline().split()[0])

data = np.zeros((num_pts,2))
for i in range(num_pts):
   (pid,x,y) = f.readline().split()
   data[i,0] = float(x)
   data[i,1] = float(y)


#np.savetxt("./data_gnuplt.dat:q:",data,fmt='%10.5f')


# compute wave frequency using dispersion relationship
def comp_freq(h,wav_num):
    grav = 9.807
    if (h<0):
        return(np.sqrt(grav*wav_num))
    else:
        return(np.sqrt(grav*wav_num*np.tanh(wav_num*h)))

def comp_eta(g,omega,pot):
    return(-1j*omega*pot/g)

# radius of cylinder, amplitude of wave, incident angle of wave
beta = 0.
rad = 1.
amp = 1.
grav = 9.807


# computed depth at z = ?,
z = 0.

#num_pts = 1

result = np.zeros((num_pts,2),dtype='complex128')


for i in range(num_pts):
    x = data[i,0]
    y = data[i,1]
    r = np.sqrt(x**2+y**2) 
    theta = np.arctan2(y,x)
    
    p1 = wk*(z+h)
    p2 = wk*h
    zfunc = np.cosh(p1)/np.cosh(p2)
    dzfunc = wk*np.sinh(p1)/np.cosh(p2)
    pot_n = 0.
    
    for m in range(31):
        eps = 2.0
        if (m == 0): eps = 1.0 

        # djn = mp.besselj(m,wk,1)
        # dyn = mp.bessely(m,wk,1)
        djn = bf.dbj(wk,m) 
        dyn = bf.dby(wk,m)
        dhn = complex(djn,dyn)        
        
        #jn =  mp.besselj(m,wk,0)
        #yn = mp.bessely(m,wk,0)

        # jnr =  mp.besselj(m,wk*rad,0)
        # ynr = mp.bessely(m,wk*rad,0)
        jnr = bf.bj(wk*rad,m)
        ynr = bf.by(wk*rad,m)
        hnr = complex(jnr,ynr)

        pot_n = pot_n-eps*complex(0,1)**m*(hnr*djn/dhn)*np.cos(m*theta)

    #----
    freq = comp_freq(h,wk)
    coef = -1j*grav*amp/freq
    # diffraction potential
    pot_dif = pot_n*coef*zfunc
    #incident wave potential
    pot_ind = coef*zfunc*np.exp(1j*wk*x)

    result[i,0] = pot_dif
    result[i,1] = pot_ind
#    print comp_eta(grav,freq,pot_dif+pot_ind)
#    print comp_eta(grav,freq,pot_dif)
#    print comp_eta(grav,freq,pot_ind)










#f = open("./tmp.txt")

#part1 = f.readlines()
#header = "TITLE = \"3D Mesh Grid Data for Element Boundary\"\n
#VARIABLES = \"X\", \"Y\"\n
#ZONE T=\"BODY MESH\" N=        481 ,E=         160 ,F=FEPOINT,ET=QUADRILATERAL



