from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from petsc4py import PETSc

from dolfinx import fem, io, mesh, plot, log, default_scalar_type, nls
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
import dolfinx
print(f"DOLFINx version: {dolfinx.__version__}")
import ufl
print(f"ufl version: {ufl.__version__}")
import matplotlib.pyplot as plt
import pyvista
from ufl import ds, dx, grad, inner
import numpy as np

##include <pybind11/stl.h>
#<pybind11/complex.h>
#<pybind11/functional.h


L=10
mesh = mesh.create_interval(comm=MPI.COMM_WORLD, points=(0,L), nx=100)
V = fem.functionspace(mesh, ("CG", 1))
gamma =float(0.1) #learning rate

#u = fem.Function(V)
#u_p = fem.Function(V)
#
#u = ufl.variable(u)
#
##Pi = dolfinx.assemble((1-u)**2)*dx
#Pi = (1-u)**2*dx
#F = ufl.diff(Pi,u)
#
#update_expr = u - gamma * dolfinx.cpp.fem.assemble_vector(F)
#
##u_p = u - gamma*dolfinx.cpp.fem.assemble_vector(F)
##
##dfn.assemble(ufl.inner(ufl.grad(trial), ufl.grad(test))) * dx

# Initialization
u = fem.Function(V)
Pi = (1-u)**2*dx

F = ufl.derivative(Pi, u)   # ufl form of the derivative

F_form = dolfinx.fem.form(F)   # Convert the ufl form to the FEniCSx weak form

F_vec = dolfinx.cpp.fem.petsc.create_vector_nest(F_form)   # Initialize a PETSc vector

 

# In the while/for loops:

with F_vec.localForm() as loc:
    loc.set(0)   # Zero out the vector components

dolfinx.fem.assemble_vector(F_vec, F_form)   # Assemble vector

#and use the following commands to update the solutions:
u_update = u.vector - F_vec   # The solution is updated and stored as a PETSc vector

 


