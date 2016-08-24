from util import util_funcs as ut
from numpy import pi
import numpy as np
from matplotlib import pyplot as pl

# x = np.linspace(0.5*pi,0.25*pi,100)
x = np.linspace(0.0*pi,2*pi,1000)
# x = np.linspace(3/4*pi-0.2,3/4*pi+0.2,1)
# x=[1/4*pi]
y1 = np.zeros_like(x)
y2 = np.zeros_like(x)

for i in range(len(x)):
    (y1[i],y2[i]) = ut.intersect_unit_square(x[i])
    
pl.scatter(y1,y2)
pl.axis('equal')
pl.grid()
pl.show()
    
