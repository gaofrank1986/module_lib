import basisfunc as bf
import numpy as np
from wave import *

max_m = 3
x = np.linspace(0,2,1000)
# y = np.zeros_like(x)
y = np.zeros([max_m+1,len(x)])
bf.var_mod.h = 10.0
bf.var_mod.wk = 2.0
bf.var_mod.s = 1.0
tmp = dispers(max_m,bf.var_mod.wk,bf.var_mod.h)
bf.var_mod.wvno1[0:max_m+1]= tmp[0:max_m+1] 
m = 0
# def f(

for j in range(max_m+1):
    for i in range(len(x)):
        # y[j,i] = bf.zm(j,x[i])
        y[j,i] = bf.ym(j,x[i])


from matplotlib import pyplot as pl


# for j in range(max_m+1):

for j in range(1,max_m+1):
    pl.plot(x,y[j,:])
pl.show()

