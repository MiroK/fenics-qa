from dolfin import *

Width, Height = 5, 10
mesh = RectangleMesh(0, 0, Width, Height, 10, 20)

class TopBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], Height)

# Your old way
markers = FacetFunctionSizet(mesh, 1)
ds = ds[markers]
TopBoundary().mark(markers, 1)
form = 1*ds(1, domain=mesh)
ans = assemble(form)   # circumnference
assert abs(ans - 2*(Width + Height)) < 1E-13
print ans

# Correct way
markers = FacetFunctionSizet(mesh, 0) # Note initialization to 0

ds = ds[markers]
TopBoundary().mark(markers, 1)
form = 1*ds(1, domain=mesh)
ans = assemble(form)   # circumnference
assert abs(ans - Width) < 1E-13
print ans
