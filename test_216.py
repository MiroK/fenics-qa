from dolfin import *

mesh = UnitSquareMesh(8, 8)
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)

def bottom(x, on_boundary): return near(x[1],0.0)
def left(x, on_boundary): return near(x[0],0.0)
bc_b = DirichletBC(W.sub(0), Constant((1.0,0.0)), bottom);
bc_l = DirichletBC(W.sub(0), Constant((0.0,0.0)), left);
bc = [bc_b, bc_l]

(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

a = inner(grad(u), grad(v))*dx + div(v)*p*dx + q*div(u)*dx
L = inner( Constant((0.0,0.0)), v)*dx

U = Function(W)
M = p * dx
solve(a == L, U, bc, tol=0.1, M=M)
