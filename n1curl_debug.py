#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
cell = tetrahedron

V = FiniteElement(cell,"N1curl", 3)

u = TestFunction(V)
v = TrialFunction(V)

s = dot(curl(u),curl(v))*dx
t = dot(u, v)*dx
