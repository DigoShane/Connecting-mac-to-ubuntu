#Here we solve the 2D Ginzbug Landau problem with an applied magnetic field.
#the whole formulation is presented in OneNote UH/superConductivity/Coding/2D Ginzburg Landau fenics.
#======================================================================================================
#The way the Code works
#1. The input to the code is:
#   a. The external field
#   b. The relaxation parameter
#   c. The absolute tolerance
#2. When reading from and writing into respective files,
#   we are writing the lagrange multiplier as a constant function
#   When reading the functions, we interpolate onto a space VAu.
#======================================================================================================
#Things to keep in mind about writing this code:-
#1. Define a functoon to evaluate the curl
#2. Define a rotation funciton.
#3. HAve replace L with l throught.
#4. All variables are lower case.

from dolfin import *
import fenics as fe
import numpy as np
from ufl import tanh
import matplotlib.pyplot as plt
import mshr


#Create mesh and define function space
l = 10
r = 0.1
domain = mshr.Rectangle(Point(-l,-l), Point(l, l)) - mshr.Circle(Point(0., 0.), r)
mesh = generate_mesh(domain1-domain2, int(l)*10)
x = SpatialCoordinate(mesh)
Va1 = FiniteElement("CG", mesh.ufl_cell(), 2)
Va2 = FiniteElement("CG", mesh.ufl_cell(), 2)
Vu = FiniteElement("CG", mesh.ufl_cell(), 2)
V = FunctionSpace(mesh, MixedElement(Va1, Va2, Vu))
Vcoord = FunctionSpace(mesh, "Lagrange", 2)#This is for ExtFile

# Define functions
da1a2u = TrialFunction(V)
(da1, da2, dr) = split(da1a2r)



#Boundary conditions
def boundary_O(x, on_boundary):
    tol = 1E-24
    return on_boundary and near(x[0], L, tol) and near(x[0],-L,tol) and near(x[1],-L,tol) and near(x[1],L,tol)
def boundary_I(x, on_boundary):
    tol = 1E-24
    return on_boundary and near(x[0]**2+x[1]**2, r, tol)

bc1 = DirichletBC(V.sub(1), 1, boundary_0) #Setting u=1 on outer boundary
bc2 = DirichletBC(V.sub(1), 0, boundary_I)
bc3 = DirichletBC(V.sub(0), 100, boundary_I)#A on the inside is 100
bcs = [bc1, bc2, bc3];


## Parameters
#kappa = Constant(1);
#Hin = input("External Magnetic field? ")
#H = Constant(Hin);
#rlx_par_in = input("relaxation parameter? ")
#rlx_par = Constant(rlx_par_in);
#tol_abs_in = input("absolute tolerance? ")
#tol_abs = Constant(tol_abs_in);
#Ae = H*x[0]
#
#
##-----------------------------------------------------------------------------------------------------------------
##!!xDx!! ##!!xDx!! Newton rhapson Approach
##-----------------------------------------------------------------------------------------------------------------
##Compute first variation of Pi (directional derivative about u in the direction of v)
##Aur = interpolate( Expression(("1","0.0", "1.5"), degree=2), V)#SC phase as initial cond.
##Aur = interpolate( Expression(("H*x[0]","0", "0"), H=H, degree=2), V)#Normal phase as initial condiiton
###Coexistence of phase as initial condition
##ul = Expression('0', degree=2, domain=mesh)
##Al = Expression('H*(0.5*Lx-x[0])', H=H, Lx=Lx, degree=2, domain=mesh)
##ur = Expression('1', degree=2, domain=mesh)
##Ar = Expression('0', degree=2, domain=mesh)
##Aur = interpolate( Expression(('x[0] <= 0.5*Lx + DOLFIN_EPS ? Al : Ar', 'x[0]<=0.5*Lx+DOLFIN_EPS ? ul : ur', '11'), ul=ul, ur=ur, Al=Al, Ar=Ar, Lx=Lx, degree=2), V)
##For 1D Vortex Solution.
##Aur = interpolate( Expression(("sqrt(2*tanh(x[0]+0.89)*tanh(x[0]+0.89)-1)","-sqrt(2)*sqrt(1-tanh(x[0]+0.89)*tanh(x[0]+0.89))","1"), degree=3), V)#1D vortex solution.
##---------------------------------------------------------------------------------------------------------------
##Reading input from a .xdmf file.
#Aur = Function(V)
#A = Function(Vcoord)
#u = Function(Vcoord)
#r = Function(RFnSp)
#data = np.loadtxt('test-2-Constraint.txt')
#y0 = data
#r = interpolate(Constant(float(y0)),RFnSp)
#A_in =  XDMFFile("test-0-Constraint.xdmf")
#A_in.read_checkpoint(A,"A",0)
#u_in =  XDMFFile("test-1-Constraint.xdmf")
#u_in.read_checkpoint(u,"u",0)
#assign(Aur,[A,u,r])
##plot(u)
##plt.title(r"$u(x)-b4$",fontsize=26)
##plt.show()
##plot(A)
##plt.title(r"$A(x)e_2-b4$",fontsize=26)
##plt.show()
#
#(A, u, r) = split(Aur)
#
#
#Scl=Constant(1.000);
#
#
#F = Scl*(-(1-u**2)*u*du + (1/kappa**2)*inner(grad(u), grad(du)) + A**2*u*du + 0.5*r*du + dr*(u-0.5) + u**2*A*dA + inner(grad(A), grad(dA)))*dx + Scl*H*dA*ds
##solver.parameters.nonzero_initial_guess = True
#solve(F == 0, Aur, bcs,
#   #solver_parameters={"newton_solver":{"convergence_criterion":"residual","relaxation_parameter":0.01,"relative_tolerance":0.000001,"absolute_tolerance":0.001,"maximum_iterations":500}})
#   solver_parameters={"newton_solver":{"convergence_criterion":"residual","relaxation_parameter":rlx_par,"relative_tolerance":0.000001,"absolute_tolerance":tol_abs,"maximum_iterations":100}})
#
#A = Aur.sub(0, deepcopy=True)
#u = Aur.sub(1, deepcopy=True)
#r = Aur.sub(2, deepcopy=True)
#
#
###Save solution in a .xdmf file
#Aur_split = Aur.split(True)
#Aur_out = XDMFFile('test-0-Constraint.xdmf')
#Aur_out.write_checkpoint(Aur_split[0], "A", 0, XDMFFile.Encoding.HDF5, False) #false means not appending to file
#Aur_out.close()
#Aur_out = XDMFFile('test-1-Constraint.xdmf')
#Aur_out.write_checkpoint(Aur_split[1], "u", 0, XDMFFile.Encoding.HDF5, False) #false means not appending to file
#Aur_out.close()
#with open("test-2-Constraint.txt", "w") as file:
#    print(float(Aur_split[2]), file=file)
#
#
#pie = assemble((1/(Lx))*((1-u**2)**2/2 + (1/kappa**2)*inner(grad(u), grad(u)) + A**2*u**2 + inner(grad(A-Ae), grad(A-Ae)))*dx )
#print("Energy density =", pie)
#Constraint = assemble( (u-0.5)*dx)
#print("Constraint violated by =", Constraint)
#
#
#plot(u)
#plt.title(r"$u(x)$",fontsize=26)
#plt.show()
#plot(A)
#plt.title(r"$A(x)e_2$",fontsize=26)
#plt.show()
#
