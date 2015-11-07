from dolfin import *

mesh = UnitSquareMesh(32, 32)

subdomain = AutoSubDomain(lambda x, on_boundary: on_boundary and near(x[0], 0))

boundary_parts = FacetFunction('size_t', mesh, 0)
subdomain.mark(boundary_parts, 2)


V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q
noslip = Constant((0, 0))


inflow = Expression(("sin(x[1]*pi/w_x)", "0.0"),w_x=3)
bc2 = DirichletBC(W.sub(0), inflow, boundary_parts, 2)

bcs = [bc2]
(v, q) = TestFunctions(W)
f = Constant((0, 0))
w = Function(W)
(u, p) = split(w)

F = (inner(grad(u), grad(v)) - div(v)*p - q*div(u))*dx-inner(f, v)*dx
a = lhs(F)
L = rhs(F)

M = inner(u, u)*dx()
tol = 1e-12
solve(a == L, w, bcs, tol = tol, M = M)
# (needed for further computation on coefficient vector)
(u, p) = w.split(True)
print "Norm of velocity coefficient vector: %.15g" % u.vector().norm("l2")
print "Norm of pressure coefficient vector: %.15g" % p.vector().norm("l2")
(u, p) = w.split()
plot(p.root_node(), title="Solution on initial mesh")
plot(p.leaf_node(), title="Solution on final mesh")
