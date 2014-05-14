from __future__ import division

from dolfin import *
import numpy as np

import matplotlib as mpl
mpl.use( "agg" )
import matplotlib.pyplot as plt
from matplotlib import animation
from math import exp
from math import sin
from math import cos
from math import pow
from math import tanh

# to shut off the output
set_log_level(PROGRESS)

nx = 20
circle = Circle (0.0,0.0, 1)

mesh = Mesh(circle, 100)
#mesh = UnitSquareMesh(nx, nx)
#mesh = UnitCircle(20)

print mesh.num_cells()


plot(mesh,axes=True)
interactive()

# definition of solving method
V = FunctionSpace(mesh, 'CG', 1)

T = 10
dt = 0.0125

# ------------------- RHS, exact solution
# RHS
def rhs(u):
  return 0

u_start = Expression("x[0] + 3")
# ------------------- RHS, exact solution

# ------------------- boundary conditions
# Define boundary segments
class outer_boundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-13   # tolerance for coordinate comparisons
        return on_boundary and  abs( x[0]*x[0] + x[1]*x[1] - 9) < tol
class left(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-13   # tolerance for coordinate comparisons
        return on_boundary and abs(x[0]) < tol


z_old = interpolate(u_start,V)

plot(z_old, interactive=True)

z = TrialFunction(V)
w = TestFunction(V)

a = inner(z,w)*dx + dt*dot(grad(z), grad(w))*dx
L = inner(z_old, w)*dx

z = Function(V)
problem = LinearVariationalProblem(a, L, z)
pdesys_Euler = LinearVariationalSolver(problem)
interactive()

# ---------- update of time dependent data
t = dt

# ---------- time steps
while t <= T+dt:

    # ---------- Newton
    # Euler
    pdesys_Euler.solve()

    # ---------- assigning of solution from previous time step
    z_old.assign(z)

    t += dt

plot(z, interactive=True, axes=True)
