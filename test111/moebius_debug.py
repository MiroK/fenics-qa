#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
moebcell = Cell('triangle', 3)
element = FiniteElement("Lagrange", cell, 1)

u = TrialFunction(element)
v = TestFunction(element)
f = Constant(triangle)

a = inner(grad(u), grad(v))*dx
L = f*v*dx
