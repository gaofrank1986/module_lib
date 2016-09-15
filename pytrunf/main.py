# from test import bj,by,dbj,dby
from spec_func import *
from numpy import *
from wave import dispers
from eigenf import eigen as e
from scipy.linalg import solve
from dutwav.vtk import ValuetoVTK

a=1.0
wk=1.0
h=10.0
bh=1.0
s=h-bh
g=9.81 # m/s^2
amp=1.0
max_wk=201
#----------------------
e.bh=bh
e.h=h
e.s=s
e.wk=wk
e.pi = pi
#===================
mmax = 5
lmax = 5
nmax = 10
# m-> 0:nmax
# j-> 0:mmax
# i-> 0:lmax

# calculate frequency(k0) using dispersion relationship
if (h<0):
    frq=sqrt(g*wk)
else:
    frq=sqrt(g*wk*tanh(wk*h))

print frq

k=zeros(max_wk+1,'float64')
wvno=dispers(max_wk,frq,h)
for i in arange(0,max_wk):
    e.wvno1[i]=wvno[i]
    k[i]=wvno[i]
k[0]=wk
''' k tesed'''

# print e.wvno1[0:max_wk+1]
# calculate zj series| eigen value for vertical eigen function
# eigen value for vertical eigen function basis Vm
def Lamda(j):
    return(j*pi/s)

for j in arange(0,mmax+1):
    e.wvzj1[j]=Lamda(j)

# print e.wvzj1[0:mmax+1]

# vertical direction eigen funcs
# def Yj(z,j):
    # if(j==0):
        # return(sqrt(2)/2)
    # else:
        # return(cos(Lamda(j)*(z+h)))

# radical direction eigen funcs
def Vm(r,j,m):
    if(j==0):
        return((r/a)**m)
    else:
        return(bessj(m,Lamda(j)*r)/bessi(m,Lamda(j)*a))

def dVm(r,j,m):
    if(j==0):
        return(m/a*(r/a)**(m-1))
    else:
        return(Lamda(j)*dbi(m,Lamda(j)*r)/bessi(m,Lamda(j)*a))


#--------Initial variable preparation finished here---------

#fmj[j,m]
FMJ=zeros((mmax+1,nmax+1),'complex128')
#amji[j,i,m]
AMJI=zeros((mmax+1,lmax+1,nmax+1),'complex128')
#ei[i,m]
emi=zeros((lmax+1,nmax+1),'complex128')
#bmij[i,j,m]
bmij=zeros((lmax+1,mmax+1,nmax+1),'complex128')

# HCCLM=DSIN(WL*S)/WL/DCOS(WL*H)/SQRT(2.0D0)
for i in arange(1,5):
    wl=k[i]
    tmp = sin(wl*s)/wl/cos(wl*h)/sqrt(2.0)
    print tmp
print "===================================="

"""
    Begin constructing coefficient matrix
"""
#FIXME : what does m=0 mean
for m in arange(0,nmax+1):

    bj0=bessj(m,k[0]*a)
    dbj0=dbj(m,k[0]*a)
    hkl0=bj0+1j*bessy(m,k[0]*a)
    dhkl0=dbj0+1j*dby(m,k[0]*a)
    '''hkl0,dhkl0 tested'''
    print "======"
    print m
    for j in arange(0,mmax+1):
        i=0
        FMJ[j,m] = 2.0*bj0/s*e.hcclm(0,j)
        '''fmj tested'''
        AMJI[j,i,m]= 2.0*hkl0/s*e.hcclm(0,j)
        for i in arange(1,lmax+1):
            bki = bessk(m,k[i]*a)
            AMJI[j,i,m] = 2.0*bki/s*e.hcclm(i,j)
            '''tested for amji '''

    i=0
    emi[0,m] = dbj0/dhkl0

    # evaluate bmij[0,0:j,m]
    for j in arange(0,mmax+1):
        bmij[0,j,m] = Lamda(j)*dVm(a,j,m)/(k[0]*dhkl0)*e.hcclm(i,j)/e.gwl(i)
        '''tested'''

    # # evaluate bmji[1:i,0:j,m]
    #FIXME sign diff
    for j in arange(0,mmax+1):
        for i in arange(1,lmax+1):
            dbk_j=dbk(m,k[j]*a)
            bmij[i,j,m] = Lamda(j)*dVm(a,j,m)/(k[i]*dbk_j)*e.hcclm(i,j)/e.gwl(i)
    '''tested'''

'''
    Finished Constructing Two coefficient part
'''
# !    [CCM] = [I]- [AMJI] [BMIJ]   

CCM=zeros((mmax+1,mmax+1,nmax+1),'complex128')
GMJ=zeros((mmax+1,nmax+1),'complex128')
bm=zeros((mmax+1,nmax+1),'complex128')
am=zeros((lmax+1,nmax+1),'complex128')
for m in arange(0,nmax+1):
    # !    [CCM] = [I]- [AMJI] [BMIJ]   
    CCM[:,:,m] = eye(mmax+1,dtype='complex128')-dot(AMJI[:,:,m],bmij[:,:,m])
    '''tested'''
    # !    [GMJ] = [FMJ]- [AMJI] [EMI]                            
    GMJ[:,m] = FMJ[:,m]-dot(AMJI[:,:,m],emi[:,m])
    '''GMJ tested'''
    # solve ax = b
    bm[:,m] = solve(CCM[:,:,m],GMJ[:,m])
    '''tested'''
    # ! --- {AM} =[BMiJ] {BM} -{EMi}
    am[:,m]=-emi[:,m]+dot(bmij[:,:,m],bm[:,m])
    '''tested,slightly different ,maybe bcoz rounding error'''
    # print m
    if m<2:
        print am[:,m]
    # print "===================================="


'''
    Compute potential
'''
from dutwav.mesh import Mesh

print "===================================="
file1="./surface_mesh.0000000.out"
m1=Mesh()
m1.read_mesh(file1,0)
pot_inc={}
pot_dif={}
t={}

for ii in m1.nodes:
    v1=zeros((nmax+1),'complex128')
    v2=zeros((nmax+1),'complex128')
    coef=zeros((nmax+1),'complex128')
    info = list(m1.nodes[ii])
    assert(info[2]<1.e-5)
    # read r
    # read z
    r=sqrt(info[0]**2+info[1]**2)
    z=info[2]
    assert(abs(info[2])<1e-5)
    if (abs(info[0]<1e-5)):
            info[0]=info[0]+1e-5*sign(info[0])

    theta=arcsin(info[1]/r)
    if(isnan(theta)):
        theta=arcsin(sign(info[1]))
    if(info[0]<0):
        theta = pi-theta

    t[ii] = theta

    for m in arange(0,nmax+1):
        eps=2.0
        if(m==0):
            eps=1.0
        bsj=bessj(m,k[0]*r)
        hkl=bsj+1j*bessy(m,k[0]*r)

        v1[m]=am[0,m]*hkl*e.zm(0,z)
        v2[m]=bsj*(e.zm(0,z))

        # # note start index
        for i in arange(1,lmax+1):
            v1[m] = v1[m]+am[i,m]*bessk(m,k[i]*r)*e.zm(m,z)

        coef[m] =-1j*g*amp/frq*eps*((1j)**m)*cos(m*theta)

    pot_dif[ii]=dot(v1,coef)
    pot_inc[ii]=dot(v2,coef)
# name1='pot_inc'+'.'+ '{:0>7d}'.format(0)
# name2='pot_dif'+'.'+ '{:0>7d}'.format(0)
# name3='theta'+'.'+ '{:0>7d}'.format(0)
# ValuetoVTK(m1,name1,pot_inc)
# ValuetoVTK(m1,name2,pot_dif)
# ValuetoVTK(m1,name3,t)
#======================================================
eta_inc={}
eta_dif={}

perd=(2*pi/frq)
maxT=3*perd
deltaT=perd/60
numIter=int(maxT/deltaT)

for jj in arange(numIter):
    print jj
    for ii in m1.nodes:
         eta_inc[ii] = (pot_inc[ii]*exp(-1j*frq*jj*deltaT)).real
         eta_dif[ii] = (pot_dif[ii]*exp(-1j*frq*jj*deltaT)).real
    name1='eta_inc'+'.'+ '{:0>7d}'.format(jj)
    name2='eta_dif'+'.'+ '{:0>7d}'.format(jj)
    ValuetoVTK(m1,name1,eta_inc)
    ValuetoVTK(m1,name2,eta_dif)
#======================================================
    





                
