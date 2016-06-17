from dolfin import *
from ufl import zero

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

cs = map(Constant, range(4))
forms = [c*inner(u, v)*dx for c in cs]

f0 = sum(forms)
f1 = sum(forms, 0)
f2 = sum(forms, zero())
f3 = sum(forms[1:], forms[0])

iforms = iter(forms)
f4 = sum(iforms, next(iforms))

print len(set([f.signature() for f in (f0, f1, f2, f3, f4)]))
