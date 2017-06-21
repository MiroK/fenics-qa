from fenics import *

mesh = IntervalMesh(10, 0, 1)

Pi = VectorElement("Lagrange",mesh.ufl_cell(),1,dim=2)
PS = MixedElement([Pi,Pi])
V = FunctionSpace(mesh,PS)

U = Function(V)
u1,u2 = split(U) 

test_u1, test_u2 = TestFunctions(V)
dU = TrialFunction(V)

foo = inner(u1,u1)
# a = conditional(gt(1,0), foo*u2[0], foo*u2[0]) #error
# a = inner(u1,u1)*u2[0] # works
a = foo*conditional(gt(1, 0), u2[0], u2[0])

b = derivative(a,u1,test_u1)

F = b*dx

problem = NonlinearVariationalProblem(F, U, [], J=derivative(F, U, dU))
