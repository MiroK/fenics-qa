import dolfin as df

# Create a unit cube mesh
mesh = df.UnitCubeMesh(10, 10, 10)

# Create a cell function to mark a region
# consisting of the left half of the mesh
region_markers = df.CellFunction('size_t', mesh)

class Domain(df.SubDomain):
    def inside(self, pt, on_boundary):
        return pt[0] <= 0.5

subdomain = Domain()
subdomain.mark(region_markers, 1)

# Create a restricted function space
restriction = df.Restriction(region_markers, 1)
V_restr = df.VectorFunctionSpace(mesh, 'CG', 1)#, restriction=restriction)
