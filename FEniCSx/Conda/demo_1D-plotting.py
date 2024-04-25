#This is a demo for 1D plotting for FEniCSx.
#This is taken from "https://fenicsproject.discourse.group/t/matplotlib-instead-of-pyvista-in-dolfinx/9733"
#---------------------------------------------------------------
# ISSUES WITH THE CODE
#1. This doesnt run. 
#   Returns a segmentation fault at plt.plot.


import dolfinx
#print(f"DOLFINx version: {dolfinx.__version__}")
import numpy as np
from mpi4py import MPI

import ufl
from dolfinx import fem, io, mesh, plot
from ufl import ds, dx, grad, inner

from petsc4py.PETSc import ScalarType


mesh = dolfinx.mesh.create_unit_interval(MPI.COMM_WORLD, 10)
P = 1
V = dolfinx.fem.FunctionSpace(mesh, ("Lagrange", P))
uh = dolfinx.fem.Function(V)
uh.interpolate(lambda x: np.sin(2*np.pi*x[0]))
x = V.tabulate_dof_coordinates()
x_order = np.argsort(x[:,0])

print(x[:,0])
print("length of x[:,0]=", len(x[:,0]))
print(x_order)
print("length of x_order=", len(x_order))
print("length of uh=", len(uh.x.array))

#print("=============================================")
#import matplotlib.pyplot as plt
#plt.plot(uh.vector.array)
##plt.plot(x[x_order,0], uh.x.array[x_order])
#plt.show()
#plt.savefig("u_1d.png")


print("=============================================")
import matplotlib.pyplot as plt
cells, types, x = plot.vtk_mesh(V)
plt.figure(1,dpi=300)  
#exact = x[:,0] - np.sinh ( x[:,0] ) / np.sinh ( 1.0 )
plt.plot ( x[:,0], uh.x.array.real)
#plt.plot ( x[:,0], exact)
filename = 'bvp_02_solution.png'
plt.savefig ( filename )
print ( '  Graphics saved as "%s"' % ( filename ) )
plt.show()
plt.close ( )


