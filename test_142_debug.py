#!/usr/bin/env python
# -*- coding: utf-8 -*-
from ufl import *
set_level(DEBUG)
cell = triangle

c = TensorConstant(triangle, (2, 2))
s = TensorConstant(triangle, (2, 2), {(1, 0):(0, 1)})

V = FiniteElement('Lagrange', cell, 1)
u = TrialFunction(V)
v = TestFunction(V)

ac = inner(dot(c, grad(u)), grad(v))*dx
as = inner(dot(s, grad(u)), grad(v))*dx
forms = [ac, as]
