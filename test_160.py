from dolfin import *
# <<<<<<< HEAD
import numpy as np

mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 2, 10)
iface = FacetFunction('size_t', mesh, 0)
CompiledSubDomain('near(x[0], 0)').mark(iface, 1)

mesh.init(1, 0)      
ivertex = list(set(sum((facet.entities(0).tolist()
                        for facet in SubsetIterator(iface, 1)), [])))

dx = np.linspace(0, 2, 10)
dx = 0.5*np.sin(2*pi*dx)

coordinates = mesh.coordinates()
for dxi in dx:
    coordinates[ivertex, 0] = dxi
    plot(mesh, interactive=True)



jvertex = []
# For each interface facet get the vertices connected to it
for facet in SubsetIterator(iface, 1):
    jvertex.extend(facet.entities(0))
# From the construction there are duplicate vertices in the list. Filter them by
# creating a set from the list
jvertex = set(jvertex)
# Bacause we want to use jvertex for indexing arrays it needs to be converted
# back to list
jvertex = list(jvertex)

print jvertex == ivertex
# =======
# mesh = Mesh('mesh_spherical.xml')    
# 
# V = FunctionSpace(mesh, "CG", 1)
# u = TrialFunction(V)
# v = TestFunction(V)
# 
# lamb = .6
# x = SpatialCoordinate(mesh)
# f = lamb*(lamb + 1.)*(sin(x[2]))**lamb*sin(lamb*x[1])
# a = Dx(u,0)*Dx(v,0)*dx
# a += 1./(x[0])**2*Dx(u,1)*Dx(v,1)*dx
# a += 1./(x[0]*sin(x[1]))**2*Dx(u,2)*Dx(v,2)*dx
# l = f*v*dx
# A = assemble(a)
# #L = assemble(l)
# >>>>>>> 556af950831cd43c4dc4cc4384d8128fcc71baa9
