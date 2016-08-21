from dutwav.mesh import Mesh

file1="./aa.TXT"
surface_type=1
maxn=999
prefix="fort.9"

m1 = Mesh()

#read mesh in external form
m1.read_mesh(file1,9)
m1._mark_surface_elems()
# define cylinder btm tag
m1._mark_elems_at_z(-1.0,'btm')
# define inner cylinder wall tag
m1._mark_elem_withr(0.98,1.1,'cwll')
# extract m3 and m2
m3 = m1.extract_mesh(['cwll'])
m3._update_tag(['extract'],'body')
m3.get_edge_info()
m3.renew_circular_midpoint()
m3._redo_bd_nrml()
m3._reverse_nrml()

m2 = m1.extract_mesh(['btm'])
m2.get_edge_info()
m2.renew_circular_midpoint()
m2._update_all_nrml([0,0,1])
m2.devour_mesh(m3)
m2.output_mesh("./inner_body.txt",0)


m4= m1.extract_mesh(['surface'])
# m4.tecplt_nrml("./model_surface.dat")

# change the midpoints to projected on the curve
m4.get_edge_info()
m4.renew_circular_midpoint()
m4.tecplt_nrml("./model_surface1.dat")

m4.output_mesh("./outer_sf.txt",1)

#FIXME this print_nodedic cannot be used becoz when outer_sf is reread ,the ordering of node will be changed 
# from dutwav.io import print_nodedic
# print_nodedic("./surface_nodes.txt",m4.nodes)

# from dutwav.util import create_animation
# create_animation("./",m1,maxn,prefix)
