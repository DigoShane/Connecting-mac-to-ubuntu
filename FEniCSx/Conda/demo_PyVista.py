from mpi4py import MPI
from petsc4py.PETSc import ScalarType  # type: ignore

import numpy as np

import dolfinx
print(f"DOLFINx version: {dolfinx.__version__}")
import ufl
print(f"ufl version: {ufl.__version__}")
from dolfinx import fem, io, mesh, plot
#from dolfinx.plot import plot
from dolfinx.fem import Function, functionspace
from dolfinx.mesh import (CellType, compute_midpoints, create_unit_cube, create_unit_square, meshtags)


try:
    import pyvista
except ModuleNotFoundError:
    print("'pyvista' is required to visualise the solution")
    exit(0)

# If environment variable PYVISTA_OFF_SCREEN is set to true save a png
# otherwise create interactive plot
if pyvista.OFF_SCREEN:
    pyvista.start_xvfb(wait=0.1)

# Set some global options for all plots
transparent = False
figsize = 800

#Plotting a finite element Function using warp by scalar
def plot_scalar():
    # We start by creating a unit square mesh and interpolating a
    # function into a degree 1 Lagrange space
    msh = create_unit_square(MPI.COMM_WORLD, 12, 12, cell_type=CellType.quadrilateral)
    V = functionspace(msh, ("Lagrange", 1))
    u = Function(V, dtype=np.float64)
    u.interpolate(lambda x: np.sin(np.pi * x[0]) * np.sin(2 * x[1] * np.pi))

    #To visualize the function u, we create a VTK-compatible grid to map
    #values of u to
    cells, types, x = plot.vtk_mesh(V)
    grid = pyvista.UnstructuredGrid(cells, types, x)
    grid.point_data["u"] = u.x.array

    #The function "u" is set as the active scalar for the mesh, and
    #warp in z-direction is set
    grid.set_active_scalars("u")
    warped = grid.warp_by_scalar()

    #A plotting window is created with two sub-plots, one of the scalar
    #values and the other of the mesh is warped by the scalar values in
    #z-direction
    subplotter = pyvista.Plotter(shape=(1, 2))
    subplotter.subplot(0, 0)
    subplotter.add_text("Scalar contour field", font_size=14, color="black", position="upper_edge")
    subplotter.add_mesh(grid, show_edges=True, show_scalar_bar=True)
    subplotter.view_xy()

    subplotter.subplot(0, 1)
    subplotter.add_text("Warped function", position="upper_edge", font_size=14, color="black")
    sargs = dict(height=0.8, width=0.1, vertical=True, position_x=0.05,
                 position_y=0.05, fmt="%1.2e", title_font_size=40, color="black", label_font_size=25)
    subplotter.set_position([-3, 2.6, 0.3])
    subplotter.set_focus([3, -1, -0.15])
    subplotter.set_viewup([0, 0, 1])
    subplotter.add_mesh(warped, show_edges=True, scalar_bar_args=sargs)
    if pyvista.OFF_SCREEN:
        subplotter.screenshot("2D_function_warp.png", transparent_background=transparent,
                              window_size=[figsize, figsize])
    else:
        subplotter.show()

































