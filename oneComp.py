from dolfin import *

name ="twocube"
subdomainMarkers=[29,229]  # from 'Physical Volume(29)' in .geo file  
boundaryMarkers = [27,227] # from 'Physical Surface(28)' in .geo file 

mesh = Mesh(name+".xml")     

## load subdomains, boundaries 
subdomains = MeshFunction("size_t", mesh, name+"_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, name+"_facet_region.xml")

dx = Measure("dx",domain=mesh)
ds = Measure("ds",domain=mesh,subdomain_data=boundaries)   

# Surfaces 
for i,marker in enumerate(boundaryMarkers):
  print "area marker ", marker, " ", assemble(Constant(1.)*ds(marker))

#print "vol", assemble(Constant(1.)*dx)
for i,marker in enumerate(subdomainMarkers):
  submesh = SubMesh(mesh, subdomains, marker)
  print "vol marker ", marker, " ", assemble(Constant(1.)*dx(submesh))


# Apply these markings to our surface/distance measures
ds = ds(domain=mesh,subdomain_data=boundaries)    
dx = dx(domain=mesh)
V = FunctionSpace(mesh,"CG",1)

v = TestFunction(V)    
u_n = Function(V)   
u_0 = Function(V) 

Du = Constant(1)
dt = 1.
#
bcs=[]
#bcs = [  
#  DirichletBC(V,Constant(2.),boundaries,boundaryMarkers[0]),
#  DirichletBC(V,Constant(0.),boundaries,boundaryMarkers[1])]


# LHS 
# works 
# F = ((u_n-u_0)*v/dt + Du*inner(grad(u_n), grad(v)))*dx()
# doesn't work 
F = ((u_n-u_0)*v/dt + Du*inner(grad(u_n), grad(v)))*dx( subdomainMarkers[0] )
F+= ((u_n-u_0)*v/dt + Du*inner(grad(u_n), grad(v)))*dx( subdomainMarkers[1] )
F+= -Constant(0.1)*v*ds( boundaryMarkers[0])
F+=  Constant(0.1)*v*ds( boundaryMarkers[1])

u_0.interpolate(Constant(2.))
u_n.vector()[:] = u_0.vector()[:]

t=0
tstop=5
while t < tstop+DOLFIN_EPS:
    solve(F==0, u_n,bcs)
    print "############### ", t

    submesh = SubMesh(mesh, subdomains, subdomainMarkers[0])
    #print "vol", assemble(Constant(1.)*dx)
    for i,marker in enumerate(subdomainMarkers):
      submesh = SubMesh(mesh, subdomains, marker)
      volSubMesh = assemble(Constant(1.)*dx(submesh))
      concSubMesh = assemble(u_n*dx(submesh))/volSubMesh
      print "conc marker ", marker, " ", concSubMesh

    u_0.assign(u_n)
    t += float(dt)


