import dolfin as df
import numpy as np

mesh = df.RectangleMesh(-1, -1, 1, 1, 2, 2)

# Create a scalar and vector spaces
space_S1 = df.FunctionSpace(mesh, 'CG', 1)
space_S3 = df.VectorFunctionSpace(mesh, 'CG', 1, dim=3)

# Create two scalar functions and a vector function
# from the previously defined spaces
f1 = df.interpolate(df.Expression("1"), space_S1)
f2 = df.Function(space_S1)
f3 = df.Function(space_S3)

# If we want a scalar function with twice
# the values of f1 (THIS WORKS)
vec = df.assemble(df.dot(2 * f1, df.TestFunction(space_S1)) * df.dP)
f2.vector().axpy(1, vec)

# Now we want to modify the vector function space
# to get the vector (2, 1, 1) in every vertex
vec = df.assemble(df.dot(df.as_vector((2 * f1, f1, f1)), df.TestFunction(space_S3)) * df.dP)
f3.vector().axpy(1, vec)

# Scalar space assign
g2 = df.Function(space_S1)
g2.assign(df.FunctionAXPY(f1, 2.))

g2.vector().axpy(-1, f2.vector())
assert df.near(g2.vector().norm('l2'), 0)

# Assigner for components of the vector space
S3_assigners = [df.FunctionAssigner(space_S3.sub(i), space_S1) for i in range(space_S3.num_sub_spaces())]

g3 = df.Function(space_S3)
# Assign to components
comps = [f2, f1, f1]
[S3_assigners[i].assign(g3.sub(i), comp) for i, comp in enumerate(comps)]

g3.vector().axpy(-1, f3.vector())
assert df.near(g3.vector().norm('l2'), 0)


df.plot(f2)
df.plot(f3)
df.interactive()
