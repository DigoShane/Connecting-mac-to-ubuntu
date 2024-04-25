from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from petsc4py import PETSc

from dolfinx.fem.petsc import create_vector, assemble_vector 

from dolfinx import fem, io, mesh, plot, log, default_scalar_type, nls
from dolfinx.fem import Function, FunctionSpace
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
import dolfinx
print(f"DOLFINx version: {dolfinx.__version__}")
import ufl
print(f"ufl version: {ufl.__version__}")
import matplotlib.pyplot as plt
import pyvista
print(pyvista.global_theme.jupyter_backend)
from ufl import ds, dx, grad, inner
import numpy as np


L=10
mesh = mesh.create_interval(comm=MPI.COMM_WORLD, points=(0,L), nx=100)
x = ufl.SpatialCoordinate(mesh)
V = FunctionSpace(mesh, ("CG", 1))
gamma = float(0.001) # learning rate
NN = int(input("Enter number of iterations: "))

# Initialization
u = Function(V)
Pi = gamma*((1-u)**2+ inner(grad(u), grad(u)) )*dx # cheat way to make sure we get \gamma F =\Delta u

F = ufl.derivative(Pi, u)   # ufl form of the derivative

F_form = dolfinx.fem.form(F)   # Convert the ufl form to the FEniCSx weak form

F_vec = create_vector(F_form)   # Initialize a PETSc vector that is compatible with a linear form

#Initialize initial conditions
u0 = fem.Function(V)
u0.interpolate( fem.Expression( x[0]/L, V.element.interpolation_points() ) ) #initial conditions
u = u0 

#In the while/for loops:
for i in range(NN):
 with F_vec.localForm() as loc:
     loc.set(0)   # Zero out the vector components. We want to use assemble_vector to initialize it below.
 
 assemble_vector(F_vec, F_form)   # Assemble vector to initialize F_form into F_vec.
 
 u_update = u.vector - F_vec   # The solution is updated and stored as a PETSc vector

 with u.vector.localForm() as loc:
   for i in range(0,u_update.size):
     loc.setValue(i,u_update[i])
 
u_out = Function(V)
u_out.vector[:] = u.vector.array


try:
    import pyvista
    cells, types, x = plot.vtk_mesh(V)
    grid = pyvista.UnstructuredGrid(cells, types, x)
    grid.point_data["u"] = u.x.array.real
    grid.set_active_scalars("u")
    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, show_edges=True)
    warped = grid.warp_by_scalar()
    plotter.add_mesh(warped)
    if pyvista.OFF_SCREEN:
        pyvista.start_xvfb(wait=0.1)
        plotter.screenshot("uh_poisson.png")
    else:
        plotter.show()
except ModuleNotFoundError:
    print("'pyvista' is required to visualise the solution")
    print("Install 'pyvista' with pip: 'python3 -m pip install pyvista'")

