#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
cell = tetrahedron
V0 = VectorFiniteElement('Lagrange', cell, 1)
V1 = VectorFiniteElement('Bubble', cell, 3)
element = V0 + V1
