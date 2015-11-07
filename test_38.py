import dolfin as dl
import numpy as np
import os

#Compile the cpp module:
cpp_sources = ["myla.cpp"]  
cpp_module = dl.compile_extension_module(
        sources=cpp_sources,
        include_dirs=["."])

#Generate a mesh and some matrices
nx = 32
ny = 32
mesh = dl.UnitSquareMesh(nx, ny)
Vh = dl.FunctionSpace(mesh, 'Lagrange', 1)

ndof = Vh.dim()

uh, vh = dl.TrialFunction(Vh), dl.TestFunction(Vh)

A = dl.assemble( dl.inner( dl.nabla_grad(uh), dl.nabla_grad(vh) ) *dl.dx )
M = dl.assemble( dl.inner( uh, vh ) *dl.dx )
ones = dl.interpolate(dl.Constant(1.), Vh).vector()
Mones = M*ones
s = Mones
s.set_local( np.ones(ndof) / Mones.array() )
M.zero()
M.set_diagonal(s)

myla = cpp_module.cpp_linalg()
B = myla.MatPtAP(M, A)

x = dl.Vector()
B.init_vector(x,1)

#This works fine
B.mult(s,x)

# Here I get an error
x.set_local(s.array())
