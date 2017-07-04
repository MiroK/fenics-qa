from dolfin import *

mesh = UnitSquareMesh(100, 10)

x = SpatialCoordinate(mesh)
f = conditional(gt(x[0], 0.5), x[0]+x[1], Constant(0))
# In the generated code f will behave as (x, y) -> 1 if x+y > 0.5 else 0
V = FunctionSpace(mesh, 'DG', 1)
fh = project(f, V)

plot(fh, interactive=True)

# Make up tensor
s = outer(grad(f), grad(f))

# Conditionals can only compare scalar valued expression
A = Constant(((0, 0), (0, 0)))
try: # Fail comparing two matrices
    conditional(eq(s, A), Constant(1), Constant(2))
except dolfin.UFLException as e:
    pass
# Compare S, A vie (S-A):(S-A)
g = conditional(le(abs(inner(s-A, s-A)), Constant(0.1)), Constant((0, 1)), Constant((0, 0)))

W = VectorFunctionSpace(mesh, 'DG', 0)
gh = project(g, W)
plot(gh, interactive=True)
