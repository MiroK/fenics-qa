from dolfin import *

element = FiniteElement("CG", tetrahedron, 1)
t = Coefficient(element)
f1 = 2*t**2
print replace(f1,{t:2}) # works

t = variable(t)
f2 = 2*t**2
df2_dt = diff(f2, t)
print replace(df2_dt, {t:2})
