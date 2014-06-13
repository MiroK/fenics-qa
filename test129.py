from dolfin import *


class NormalVector(Expression):
    def __init__(self, mesh, facet_f):
        self.mesh = mesh
        self.facet_f = facet_f
        tdim = mesh.topology().dim()
        self.mesh.init(tdim, tdim-1)
        self.bb_tree = mesh.bounding_box_tree()
        self.bb_tree.build(mesh, tdim-1)
        self.tdim = tdim

    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        facets = cell.entities(self.tdim-1)
        bdr_facets = [facet
                      for facet in facets if Facet(self.mesh, facet).exterior()]
        print 'cell index', cell.index()
        print 'bdr facets', bdr_facets
        print 'x', x
        for bdr_facet in bdr_facets:
            print 'isect', self.bb_tree.compute_entity_collisions(Point(*x))
        print
        values[0] = 0
        values[1] = 0

    def value_shape(self):
        return (2, )

mesh = UnitSquareMesh(2, 2) #CircleMesh(Point(0., 0.), 1, 0.5)
facet_f = FacetFunction('bool', mesh, False)
DomainBoundary().mark(facet_f, True)

V = VectorFunctionSpace(mesh, 'CG', 1)
n = NormalVector(mesh, facet_f)
interpolate(n, V)
