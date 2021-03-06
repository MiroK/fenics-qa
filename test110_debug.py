#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
#kind of element
gdim = 2
cell = tetrahedron if gdim == 3 else triangle

# Function Space for Velocity U and Pressure V 
U = VectorElement("CG", cell, 2)
V = FiniteElement("CG", cell, 1)

# Mixed Function Space 
W = U*V

# (initial guess) Function for solution
w = Coefficient(W)

# (right hand side) Function
f = Coefficient(U);

# split into parts, w_u is desired solution
(w_u,w_p) = (as_vector(tuple([w[i] for i in range(gdim)])), w[gdim])

# TestFunction
(psi_u,psi_p) = TestFunctions(W)

# user defined viscosity
# NU = Constant(cell)
# Weak Formulation
# nonlinear termin (+) 
# laplace (+)
# pressure (+)  
# solenoidal (+) 
# rhs
 F =  inner(psi_u,dot(grad(w_u),w_u))*dx + \
     inner(grad(psi_u),grad(w_u))*dx    + \
     -1*div(psi_u)*w_p*dx               + \
     -1*psi_p*div(w_u)*dx               + \
     dot(psi_u,f)*dx                          

 # Derivative
 dw = TrialFunction(W)
 J = derivative(F, w, dw)
