from dolfin import *

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)

u, v = TrialFunction(V), TestFunction(V)
k = Constant(1)
dx = dx(domain=mesh)

try:
    L = v*u*(dx(1)+(1-k)*(dx(2)+dx(3)))
except TypeError as e:
    print e
    first = dx(1)
    print 'dx(1) is', type(first)
    second = (1-k)*(dx(2)+dx(3))
    print '(1-k)*(dx(2)+dx(3))', type(second)
    # You could try to change first so that it becomes form, e.g.
    first = 1*dx(1)
    bracket = first + second
    # Howver
    try:
        L = v*u*bracket
    except TypeError as e:
        print e
        # And so
        L = v*u*dx(1) + v*u*(1-k)*(dx(2)+dx(3))
