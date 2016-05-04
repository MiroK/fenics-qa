#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
cell = triangle

S = TensorElement('Lagrange', cell, 1, symmetry=True)
W = TensorElement('Lagrange', cell, 1, symmetry=False)
elements = [S, W]
