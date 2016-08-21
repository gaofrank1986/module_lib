from dutwav.mesh import Mesh

file1="./outer_sf.txt"
# file2="./INPUT/DATWPMS.txt"
surface_type=1
maxn=999
prefix="fort.9"

m1 = Mesh()
m2 = Mesh()

# m1.read_mesh(file1,surface_type)
# m1.read_mesh(file1,0)
# m1.draw_model()
from dutwav.vtk import MeshtoVTK
m1.read_mesh("./inner_body.txt",0)
MeshtoVTK(m1,'cylinder_btm')
for i in m1.nodes:
    info=list(m1.nodes[i])
    info[2]=-3.*info[2]
    m1.nodes[i]=tuple(info)
MeshtoVTK(m1,'cylinder')


m2.read_mesh(file1,surface_type)
m2.read_mesh(file1,1)

# Get surface node info
from dutwav.io import print_nodedic
print_nodedic("./nodedic.txt",m2.nodes)

# m2.tecplt_nrml("./model_nrml.dat")

# from dutwav.util import create_animation
# create_animation("./",m1,maxn,prefix)
from dutwav.postprocess import display_value
# display_value("./fort.100","./inc2_pot.dat",m2)

# display_value("./test_solid_angle.txt","./solid_angle.dat",m2)
# display_value("./out_sla.txt","./solid_angle.dat",m2)


### read file result to output in vtk

# from dutwav.vtk import waveVtk
# waveVtk(m2,"test","./fort.100")

from dutwav.result import AnalyticalPotential
# a=AnalyticalPotential()
# a.read_vector("./fort.99")
# a.w=0.5
# a.generate_wave(m2,24,24./200,factor=6.,prefix='total_wav')

a=AnalyticalPotential()
a.read_vector("./fort.100")
a.w=0.5
a.generate_wave(m2,24,24./200,factor=6.,prefix='inc_wav')

a=AnalyticalPotential()
a.read_vector("./fort.101")
a.w=0.5
a.generate_wave(m2,24,24./200,factor=6.,prefix='dif_wav')
