from dutwav.mesh import Mesh

file1="./inner_top.txt"
surface_type=1
maxn=999
prefix="fort.9"

m1 = Mesh()
m2 = Mesh()

m1.read_mesh(file1,0)
# m1.read_mesh(file1,0)
# m1.draw_model()


m1.tecplt_nrml("./model_nrml.dat")
# m1.tecplt_nrml("./model_nrml2.dat",3)


# from dutwav.util import create_animation
# create_animation("./",m1,maxn,prefix)
