import bessel_func as bsf
import numpy as np

max_m = 3
x = np.linspace(0.0001,10,1000)
# y = np.zeros_like(x)
y = np.zeros([max_m+1,len(x)])
# bf.var_mod.h = 10.0
# bf.var_mod.wk = 2.0
# bf.var_mod.s = 1.0
# tmp = dispers(max_m,bf.var_mod.wk,bf.var_mod.h)
# bf.var_mod.wvno1[0:max_m+1]= tmp[0:max_m+1] 


for j in range(max_m+1):
    print j
    for i in range(len(x)):
        y[j,i] = bsf.bk(x[i],j)
        # y[j,i] = bsf.by(x[i],j)
        # y[j,i] = bsf.bk(x[i],j)
        # y[j,i] = bsf.bi(x[i],j)


from matplotlib import pyplot as pl

ls =['o-','.-','+-','x-']
# for j in range(max_m+1):

for j in range(max_m+1):
    pl.plot(x,y[j,:])
    # pl.plot(x,y[j,:],ls[j])
# pl.axis([0,10,-2,1])
pl.axis([0,5,0,2])
pl.grid()
# pl.legend((line1,line2,line3,line4),'m=0','m=1','m=2','m=3')
pl.legend(['m=0','m=1','m=2','m=3'])
pl.show()



