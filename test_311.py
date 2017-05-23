from dolfin import *
import numpy as np

mesh = UnitCubeMesh(10, 2, 3)
V = FunctionSpace(mesh, 'N1curl', 1)  

# Let's conjecture that given edge (v0_index, v1_index) fenics sets the tangent
# in direction of vertices[v1_index]-vertices[v0_index]. Now for N1curl the dofs
# are int_edge edge_tangent*u. If we pass our tangent as u we shall see by the sign
# of dof evaluation if the tangent is right

Vdmap = V.dofmap()
x = mesh.coordinates()

u = Expression(('t0', 't1', 't2'), t0=0, t1=0, t2=0, degree=0)

mesh.init(3, 1)
mesh.init(1, 0)
elm = V.dolfin_element()
for cell in cells(mesh):
    cell_dofs = Vdmap.cell_dofs(cell.index())
    for ei, edge in enumerate(edges(cell)):
        v0, v1 = edge.entities(0)
        tau = x[v1]-x[v0]
        tau /= np.linalg.norm(tau)

        u.t0, u.t1, u.t2 = tau[0], tau[1], tau[2]

        sign = (1./edge.length())*elm.evaluate_dof(ei, 
                                                   u,
                                                   cell.get_vertex_coordinates(),
                                                   cell.orientation(),
                                                   cell)
        assert near(sign, 1.0, 1E-12)
