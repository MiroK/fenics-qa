#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
cell = tetrahedron
element = VectorElement("Lagrange", cell, 1 )
property = FiniteElement("DG", cell, 0 )

v = TestFunction(element)
u = TrialFunction(element) 
f = Coefficient(element)

E = Coefficient(property)
nu = Coefficient(property)
epsr = Coefficient(property)

mu = E / (2*(1+nu))
lmbda = E*nu / ((1+nu)*(1-2*nu))

def epsilon(v):
    return 0.5*(grad(v) + transpose(grad(v)))

def sigma(v):
    return 2*mu*epsilon(v) + lmbda*(tr(epsilon(v)) * Identity(cell.d))

g = (2*mu + 3*lambda)*epsr*Identity(cell.d)

a = inner( sigma(u), grad(v) )*dx
L = dot( f, v ) * dx + inner(g, grad(v))*dx
