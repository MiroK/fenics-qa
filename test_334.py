from dolfin import *
import os
import numpy as np

mesh = RectangleMesh(Point(0.0, 0.0), Point(1.0, 4.0), 10, 40,"crossed")

center = Point(0.5, 0.5)
radius = 0.2

class InitialCondition(Expression):
    def eval(self, values, x):
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0
        values[3] = sqrt((x[0]-center[0])*(x[0]-center[0]) + (x[1]-center[1])*(x[1]-center[1]))-radius
    def value_shape(self):
        return (4,)

V=VectorElement("Lagrange", mesh.ufl_cell(), 2)
P=FiniteElement("Lagrange", mesh.ufl_cell(), 1)
L=FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W=MixedElement([V,P,L])

S0 = FunctionSpace(mesh,V)
S1 = FunctionSpace(mesh,P)
S2 = FunctionSpace(mesh,L)
W0 = FunctionSpace(mesh,W)

ic = InitialCondition(element=W0.ufl_element())
w = interpolate(ic, W0)
