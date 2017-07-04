from fenics import *

# define auxiliary functions P^{lm}_k(y) = y_m*delta[k,l]
class Auxiliaryfunc(Expression):
      def __init__(self, l,  m, **kwargs):
         self.l = l     
         self.m = m
         return None
      def eval(self, values, x):
         if self.l == 1:
            values[0] = x[self.m-1] 
            values[1] = 0.0
         elif self.l == 2:
            values[0] = 0.0 
            values[1] = x[self.m-1]
         else:
            values[0] = 1000.0  
            values[1] = 1000.0
      def value_shape(self):
         return (2,)

# define tensor
def eps(v):
    return as_vector([v[0].dx(0), 
                      v[1].dx(1),
                      (v[0].dx(1) + v[1].dx(0))])

mesh = RectangleMesh(Point(0.0,0.0),Point(1.0,1.0),30,30,"left/right")
V  = VectorFunctionSpace(mesh, "CG", 1)

# define bc
def boundary(x, on_boundary):
    return on_boundary

bcs = DirichletBC(V, Constant((0.0, 0.0)), boundary)

# define functions 
p = Auxiliaryfunc(1, 1, degree=1, #element=V.ufl_element())
p = interpolate(p, V)
v = TestFunction(V)
u = TrialFunction(V) 

# define and solve variational problem
C = as_matrix([[5.0, 1.0, 0.],[1.0, 4.0, 0.],[0., 0., 2.0]]) 
a = inner((C * eps(u)), eps(v))*dx
L= inner((C * eps(p)), eps(v))*dx
assemble(L)

w = Function(V)    
Problem = LinearVariationalProblem(a, L, w, bcs)
Solver  = LinearVariationalSolver(Problem)

Solver.solve()

plot(w, interactive=True)
