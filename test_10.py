#!/usr/bin/env python
#
#  Use finite elements to solve Poisson equation 
#  with Robin and periodic boundary conditions 
#
#


from dolfin import *           
from mshr import *        


#  Define domain

side     =  2.
radius   =  0.35
res      =  20                   #  number of nodes along boundary side     

radius2  =  radius*radius
halfside =  0.5 * side
tol      =  side/res/100.        #  tolerance for coordinate comparisons

domain   =  Rectangle(dolfin.Point(-halfside, -halfside),  
                      dolfin.Point(halfside,halfside))     -   \
            Circle(dolfin.Point(0.0, 0.0), .35)


#  Define source and fluxes

q = Constant(  0)                #  for the moment, no production on the border
p = Constant( 10)                #  degradation on the border
f = Constant( .1)                # production in the bulk


#  Define boundaries

class InnerCircle(SubDomain):
   def inside(self, x, on_boundary):
      return on_boundary and (x[0]*x[0]+x[1]*x[1] - radius2 < tol)


rbc = InnerCircle()


class SquareBorder(SubDomain):

   def inside(self, x, on_boundary):
      # return True if on left or bottom boundary 
      # AND NOT on one of the two corners (-L/2, L/2) and (L/2, -L/2)
      return (on_boundary and 
              ( (abs(x[0]+halfside)<tol) or  (abs(x[1]+halfside)<tol) ) and not
            ( ( (abs(x[0]+halfside)<tol) and (abs(x[1]-halfside)<tol) ) or 
              ( (abs(x[0]-halfside)<tol) and (abs(x[1]+halfside)<tol) )) ) 

   def map(self, x, xn):
      # identify points on the opposite sides
      if ( (abs(x[0]-halfside)<tol) and (abs(x[1]-halfside)<tol) ):
         xn[0] = x[0] - side
         xn[1] = x[1] - side
      elif ( abs(x[0]-halfside)<tol ):
         xn[0] = x[0] - side
         xn[1] = x[1]
      else:   # abs(x[1]-halfside)<tol )
         xn[0] = x[0]
         xn[1] = x[1] - side


pbc = SquareBorder()


#  Create mesh and define function space



mesh=generate_mesh(domain,res)


foo = PeriodicBoundaryComputation().masters_slaves(mesh, pbc, 1)
plot(foo, interactive=True)

print 'befire'
V = FunctionSpace(mesh, "Lagrange", 1, constrained_domain=pbc)
print 'after'


left = AutoSubDomain(lambda x, on_boundary: on_boundary and near(x[0], -halfside))
right = AutoSubDomain(lambda x, on_boundary: on_boundary and near(x[0], halfside))
top = AutoSubDomain(lambda x, on_boundary: on_boundary and near(x[1], halfside))
bottom = AutoSubDomain(lambda x, on_boundary: on_boundary and near(x[1], -halfside))

facet_f = FacetFunction('size_t', mesh, 0)
left.mark(facet_f, 1)
right.mark(facet_f, 2)
top.mark(facet_f, 3)
bottom.mark(facet_f, 4)

left_facets = set(f.index() for f in SubsetIterator(facet_f, 1))
right_facets = set(f.index() for f in SubsetIterator(facet_f, 2))
top_facets = set(f.index() for f in SubsetIterator(facet_f, 3))
bottom_facets = set(f.index() for f in SubsetIterator(facet_f, 4))

mesh.init(0, 1)
left_heights = []
right_heights = []
top_heights = []
bottom_heights = []
for vertex in vertices(mesh):
    if len(set(vertex.entities(1)) & left_facets):
        left_heights.append(vertex.x(1))

    if len(set(vertex.entities(1)) & right_facets):
        right_heights.append(vertex.x(1))

    if len(set(vertex.entities(1)) & top_facets):
        top_heights.append(vertex.x(0))

    if len(set(vertex.entities(1)) & bottom_facets):
        bottom_heights.append(vertex.x(0))

left_heights = sorted(left_heights)
right_heights = sorted(right_heights)
top_heights = sorted(top_heights)
bottom_heights = sorted(bottom_heights)

print left_heights
print right_heights
print top_heights
print bottom_heights




# Define boundary segments for Robin and periodic boundary conditions

boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

rbc.mark(boundary_parts, 0)
pbc.mark(boundary_parts, 1)


# Define variational problem

u = TrialFunction(V)
v = TestFunction(V)
dss=ds(domain=mesh, subdomain_data=boundary_parts)
a = inner(nabla_grad(u), nabla_grad(v))*dx + p*u*v*dss(0)
l = f*v*dx + p*q*v*dss(0)           


# Assemble conditions

A = assemble(a)
L = assemble(l)


# Solve A U = b 
#
#  lu    :  LU factorization (Gaussian elimination)
#  ilu   :  incomplete LU factorization (preconditioner)
#  cg    :  conjugate gradient (Krylov solver)

set_log_level(PROGRESS)
u = Function(V)
U = u.vector()
#solve(A,U,L)    
solve(A,U,L,'lu')    
#solve(A,U,L,'cg','ilu') 


#  Compute nodal values

u_ar = u.vector().array()
coor = mesh.coordinates()


# Plot solution and mesh
plot(u)
offset=-u_ar.min()+0.5*(u_ar.max()-u_ar.min())
print offset
plot(u+offset)
interactive()
