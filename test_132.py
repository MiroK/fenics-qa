from dolfin import *
parameters["form_compiler"]["representation"] = "uflacs"

mesh = UnitSquareMesh(2,2)
VV = VectorFunctionSpace(mesh, "CG", 1)
V = VV.sub(0).collapse()

r2 = 1.0/sqrt(2)
e = interpolate(Constant([r2,r2]), VV)

phi = project(atan_2(e[1], e[0]), V)
plot(phi, interactive = True)
