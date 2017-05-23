from __future__ import print_function
from dolfin import *

comm    = mpi_comm_world()
mpiRank = MPI.rank(comm)
parameters["ghost_mode"] = "shared_facet"

mesh = UnitSquareMesh(10,10)
DG0_Element = FiniteElement('DG',mesh.ufl_cell(),0)
V = FunctionSpace(mesh,DG0_Element)
u = TrialFunction(V)
v = TestFunction(V)

cell_f = MeshFunction('size_t', mesh, 0)
if mpiRank == 0:
    cell_f[0] = 1
dX = Measure('dx', domain=mesh, subdomain_data=cell_f, subdomain_id=1)

# just an example
a = inner(jump(u),jump(v))*dS + inner(u, v)*dX
L = Constant(1.0)*v*dx

# Question: how to construct a boundary condition in order to remove
# the non-trivial null-space from the problem? 
bcs = []

A,b = assemble_system(a,L,bcs)

null_space = project(Constant(1.0),V)
# test the null-space
tmp = Function(V) 
A.mult(null_space.vector(),tmp.vector())
tmp_norm = norm(tmp.vector(),'l2')
if mpiRank==0:
    print('norm = %e' % tmp_norm)
