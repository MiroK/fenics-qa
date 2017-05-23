from dolfin import *

class FacetNormalExpr(Expression):
    def __init__(self, bdry_facets, **kwargs):
        self.bdry_facets = bdry_facets

    def eval_cell(self, values, x, cell):
        mesh = self.bdry_facets.mesh()
        cell = Cell(mesh, cell.index)

        for facet in facets(cell):
            values[:] = 0
            if self.bdry_facets[facet]:
                n = facet.normal()
                values[0] = n[0]
                values[1] = n[1]
                break

    def value_shape(self): return (2, )


mesh = UnitSquareMesh(20, 20)
facet_f = FacetFunction('bool', mesh, False)
DomainBoundary().mark(facet_f, True)

n = FacetNormalExpr(bdry_facets=facet_f, degree=1)
plot(n, mesh=mesh)
interactive()

