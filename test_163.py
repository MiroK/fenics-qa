from dolfin import *

mesh = UnitSquareMesh(10, 10)

f = Constant(10)
v = f*ds(domain=mesh) + (f**2)*ds(domain=mesh) + (f**3)*ds(domain=mesh)

for integral in v.integrals():
    integrand = integral.integrand()
    form = integrand*dx(domain=mesh)

    print form, assemble(form)
