#!/usr/bin/env python
# -*- coding: utf-8 -*-
from ufl import *
set_level(DEBUG)
cell = triangle

V = VectorElement("Lagrange", cell, 2)
v = TestFunction(V)

W = FiniteElement("DG", cell, 0)
p = Coefficient(W)

# Does not work with Illegal component error
b = conditional(lt(1,p), VectorConstant(cell), 2*VectorConstant(cell))

# Fix
b0 = conditional(lt(1, p), Constant(cell), 2*Constant(cell))
b1 = conditional(lt(1, p), Constant(cell), 2*Constant(cell))
b = as_vector((b0, b1))

def vector_conditional(predicate, tvalue, fvalue):
    assert tvalue.ufl_shape == fvalue.ufl_shape
    
    cond = []
    for i in range(tvalue.ufl_shape[0]):
        ans.append(conditional(predicate, tvalue[i], fvalue[i]))
    return as_vector(cond)

b = vector_conditional(lt(1, p), VectorConstant(cell), 2*VectorConstant(cell))

L = inner(b, v)*dx
