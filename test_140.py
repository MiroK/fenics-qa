from dolfin import *
parameters["ghost_mode"] = "shared_facet"

comm = mpi_comm_world()
mpirank = MPI.rank(comm)
mpisize = MPI.size(comm)

ny = 5
mesh = UnitSquareMesh(ny-1, ny-1)
meshtopo = mesh.topology()
V = FunctionSpace(mesh, 'Lagrange', 1)

interfaces = FacetFunction('size_t', mesh)
interfaces.set_all(0)

shared_edges = meshtopo.shared_entities(1)
for edge in shared_edges.keys():
    interfaces[edge] = 1

plot(interfaces)
interactive()

dS = Measure('dS', domain=mesh, subdomain_data=interfaces) 

u = TrialFunction(V)
v = TestFunction(V)

mrobin =  inner(u('-'),v('-'))*dS(1)
Mrobin = assemble(mrobin)
print Mrobin.array()
