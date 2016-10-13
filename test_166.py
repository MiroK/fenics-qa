from dolfin import *

mesh = UnitSquareMesh(10, 10)

V = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, V)
W = TensorElement('Lagrange', mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, W)

vs = [interpolate(Constant(i+1), V) for i in range(4)]
w = Function(W)
assign(w, vs)

for i in range(2):
    for j in range(2):
        print sqrt(assemble(inner(w[i, j], w[i, j])*dx))
